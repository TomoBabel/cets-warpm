"""Warp tilt series -> CETS ``Region`` (+ companion entries)."""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
from cryoet_alignment.io.cets import ctf as cets_ctf
from cryoet_alignment.io.cets.alignment import ReferenceVolume, alignment_to_cets
from cryoet_alignment.io.cets.cli_support import Gate, SeriesReport
from cryoet_alignment.io.cets.companion import (
    AlignmentCompanion,
    ImageCompanion,
    TiltSeriesCompanion,
    TomogramCompanion,
)
from cryoet_alignment.io.cets.config import Resolver
from cryoet_alignment.io.cets.entities import (
    movie_stack_series_entity,
    region_entity,
    tilt_series_entity,
    tomogram_entity,
)
from cryoet_alignment.io.cets.frames import FRAME_CONVENTIONS, image_frame
from cryoet_alignment.io.cryoet_data_portal import Alignment
from cryoet_alignment.io.warp import WarpAlignment

from cets_warpm.discover import WarpSeries, mrc_header

ALIGNMENT_NAME = "warp"


@dataclass
class SeriesResult:
    region: Any
    tilt_series_companion: TiltSeriesCompanion
    alignment_companion: AlignmentCompanion
    tomogram_companions: Dict[str, TomogramCompanion]


def _rel(path: Optional[Path], base: Path, mode: str) -> Optional[str]:
    if path is None:
        return None
    return str(Path(path).resolve()) if mode == "absolute" else os.path.relpath(Path(path).resolve(), base.resolve())


def warp_to_cets(s: WarpSeries, res: Resolver, sr: SeriesReport, *, out_dir: Path) -> SeriesResult:
    stem = s.stem
    pix = res.require("pix", discovered=s.get("pix"), note=s.source("pix"))
    image_dims_a = res.require("image_dims_a", discovered=s.get("image_dims_a"), note=s.source("image_dims_a"))
    volume_dims_a = res.require("volume_dims_a", discovered=s.get("volume_dims_a"), note=s.source("volume_dims_a"))
    paths_mode = res.value("paths", default="relative")
    drop_locals = bool(res.optional("drop_locals", absent=False))
    no_ctf = bool(res.optional("no_ctf", absent=False))

    w = WarpAlignment.from_file(s.xml_path, pixel_size_a=pix, image_dims_a=list(image_dims_a), volume_dims_a=list(volume_dims_a))
    if w.grid_audit.has_varying_grids and not drop_locals:
        raise ValueError(
            f"{stem}: spatially varying deformation grids ({', '.join(w.grid_audit.varying_grid_names)}); "
            "the rigid profile cannot represent them - pass --drop-locals to keep the rigid part only",
        )
    if w.grid_audit.has_varying_grids:
        sr.dropped.append(f"local deformation grids {w.grid_audit.varying_grid_names}")
    if w.grid_audit.has_angle_grids:
        sr.dropped.append(f"non-zero {w.grid_audit.angle_grid_names} (affect particle orientations; not modelled)")
    if any(abs(e.movement_x) > 0 or abs(e.movement_y) > 0 for e in w.entries) or any(w.grid_audit.constant_volume_warp):
        sr.gates.append(Gate("constant_grids_folded", True,
                             value={"movement_max_a": max(max(abs(e.movement_x), abs(e.movement_y)) for e in w.entries),
                                    "volume_warp_a": w.grid_audit.constant_volume_warp}))
    hub = Alignment.from_warp(w, allow_varying_grids=drop_locals)

    n = w.n_tilts
    width, height = (int(round(image_dims_a[0] / pix)), int(round(image_dims_a[1] / pix)))
    if s.get("image_dims_px"):
        width, height = s.get("image_dims_px")
    # Warp <Angles> are -(stage angle); the nominal stage angle is their negative
    nominal = [-e.tilt_angle for e in w.entries]
    doses = [e.dose for e in w.entries]  # accumulated before the image (exclusive), Warp convention
    # acquisition index = rank of the accumulated dose (only meaningful when the doses are all distinct)
    order = [int(r) + 1 for r in np.argsort(np.argsort(doses, kind="stable"), kind="stable")] if len(set(doses)) == n else None
    if order is None:
        sr.warnings.append("<Dose> values are not distinct: acquisition order not derivable")
    exposure_const = res.optional("dose_per_tilt", discovered=s.get("exposure_per_tilt"), note=s.source("exposure_per_tilt"))

    ctfs = None
    if w.has_ctf and not no_ctf:
        ctfs = [cets_ctf.from_warp_values(e.defocus_um, e.defocus_delta_um, e.defocus_angle_deg, e.phase_shift_pi) for e in w.entries]
        if w.grid_audit.ctf_grid_spatial:
            sr.warnings.append("CTF grids vary spatially; per-tilt means used")
    movie_paths = [e.movie_path for e in w.entries]
    have_movies = all(movie_paths)
    movie_ids = [f"{stem}_movie_{i}" for i in range(n)] if have_movies else None
    image_paths = None
    if have_movies and s.tomostar_dir:
        image_paths = []
        for mp in movie_paths:
            movie = (s.tomostar_dir / mp.replace("\\", "/")).resolve()
            avg = movie.parent / "average" / (movie.stem + ".mrc")
            image_paths.append(_rel(avg, out_dir, paths_mode) if avg.exists() else None)
        if any(p is None for p in image_paths):
            image_paths = None
            sr.warnings.append("frame averages (average/<name>.mrc) not found for every tilt: TiltImage.path left null")

    ts = tilt_series_entity(
        tilt_series_id=stem, path=None, width=width, height=height, pixel_size_a=pix, nominal_angles=nominal,
        doses=doses, ctfs=ctfs, movie_stack_ids=movie_ids, image_paths=image_paths,
        movie_stack_series_id=f"{stem}_movies" if have_movies else None,
    )
    movie_series = []
    if have_movies:
        stacks = []
        for i, mp in enumerate(movie_paths):
            movie = (s.tomostar_dir / mp.replace("\\", "/")).resolve() if s.tomostar_dir else Path(mp)
            stacks.append({"id": movie_ids[i], "path": _rel(movie, out_dir, paths_mode)})
        movie_series.append(movie_stack_series_entity(series_id=f"{stem}_movies", stacks=stacks))

    # reference volume: Warp's native box at the tilt-image pixel (bin-1 grid of the reconstruction frame)
    vol_px = tuple(int(round(v / pix)) for v in volume_dims_a)
    ref_tomo = tomogram_entity(tomogram_id=f"{stem}_volume", path=None, size_px=vol_px, voxel_size_a=pix, tilt_series_id=stem)
    tomograms = [ref_tomo]
    tomo_comp = {ref_tomo.id: TomogramCompanion(voxel_implied_a=pix, source_ref="Warp VolumeDimensionsAngstrom / Tomo.Dimensions box; no file")}
    for rec in s.reconstructions:
        try:
            h = mrc_header(rec)
        except Exception as e:  # noqa: BLE001
            sr.warnings.append(f"{rec.name}: unreadable header ({e})")
            continue
        implied = volume_dims_a[0] / h["nx"]
        tomo = tomogram_entity(tomogram_id=f"{stem}_tomo_{h['voxel'][0]:.3f}", path=_rel(rec, out_dir, paths_mode),
                               size_px=(h["nx"], h["ny"], h["nz"]), voxel_size_a=implied, tilt_series_id=stem)
        tomograms.append(tomo)
        tomo_comp[tomo.id] = TomogramCompanion(voxel_header_a=h["voxel"][0], voxel_implied_a=implied,
                                               reconstruction_software="Warp", source_ref=rec.name)

    cets_alignment = alignment_to_cets(
        hub, tilt_series_id=stem, alignment_name=ALIGNMENT_NAME, image=image_frame(ts.images[0]),
        reference=ReferenceVolume.from_tomogram(ref_tomo), frame=FRAME_CONVENTIONS["WARP"],
    )
    kept = [e.z_index for e in w.entries if e.use_tilt]
    sr.gates.append(Gate("rows", len(cets_alignment.projection_alignments) == len(kept), value=len(cets_alignment.projection_alignments), expected=len(kept)))

    images = {}
    for i, e in enumerate(w.entries):
        images[ts.images[i].id] = ImageCompanion(
            acquisition_index_1b=None if order is None else order[i],
            exposure_dose=None if exposure_const is None else float(exposure_const),
            stage_angle_deg=nominal[i],
            frame_name=Path(e.movie_path.replace("\\", "/")).name if e.movie_path else None,
            use_tilt=e.use_tilt,
        )
    ts_comp = TiltSeriesCompanion(
        source_tool="Warp",
        voltage_kv=res.optional("voltage", discovered=s.get("voltage"), note=s.source("voltage")),
        cs_mm=res.optional("cs", discovered=s.get("cs"), note=s.source("cs")),
        amplitude_contrast=res.optional("amp_contrast", discovered=s.get("amp_contrast"), note=s.source("amp_contrast")),
        dose_rate=None if exposure_const is None else float(exposure_const),
        pixel_size_acquisition_a=s.get("pixel_size_acquisition_a"),
        pixel_size_ctf_a=w.ctf_pixel_size_a,
        tilt_axis_nominal_deg=float(np.median([e.tilt_axis_angle for e in w.entries if e.use_tilt])),
        alpha_offset_deg=-w.level_angle_y,
        are_angles_inverted=w.are_angles_inverted,
        images=images,
    )
    dropped = list(sr.dropped)
    aln_comp = AlignmentCompanion(
        name=ALIGNMENT_NAME, tilt_series_id=stem, format="WARP", alignment_type="GLOBAL", is_portal_standard=True,
        reference_tomogram_id=ref_tomo.id, tomogram_ids=[t.id for t in tomograms],
        native_volume_dimension_a=hub.volume_dimension, frame_convention={"image_center": "half", "volume_center": "half"},
        dropped=dropped, source_ref=s.xml_path.name,
    )
    region = region_entity(region_id=stem, tilt_series=[ts], alignments=[cets_alignment], tomograms=tomograms,
                           movie_stack_series=movie_series)
    sr.provenance = res.provenance()
    return SeriesResult(region, ts_comp, aln_comp, tomo_comp)
