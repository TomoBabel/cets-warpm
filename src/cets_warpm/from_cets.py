"""CETS region -> Warp project (``warp_tiltseries.settings``, ``tomostar/<stem>.tomostar``,
``warp_tiltseries/<stem>.xml``). Frames are not written; the hint says where Warp expects them."""

import os
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
from cryoet_alignment.io.cets import ctf as cets_ctf
from cryoet_alignment.io.cets.alignment import (
    ReferenceVolume,
    alignment_from_cets,
    alignment_name_of,
    fold_projection,
    select_alignment,
    select_tomogram,
)
from cryoet_alignment.io.cets.cli_support import Gate, SeriesReport
from cryoet_alignment.io.cets.companion import Companion
from cryoet_alignment.io.cets.config import Resolver
from cryoet_alignment.io.cets.frames import FRAME_CONVENTIONS, find_by_id, image_frame
from cryoet_alignment.io.cryoet_data_portal import Alignment
from cryoet_alignment.io.warp import TomostarRow, WarpAlignment, WarpSettings, WarpTomostar
from cryoet_alignment.io.warp.alignment import _CTF_GRIDS  # noqa: F401 - documents the grids written

SETTINGS = "warp_tiltseries.settings"
TOMOSTAR_DIR = "tomostar"
PROCESSING_DIR = "warp_tiltseries"
OPERATOR_TOL_DEG = 1e-6


def _conflicts(existing: WarpSettings, wanted: WarpSettings) -> List[str]:
    out = []
    for section, name in (("Import", "PixelSize"), ("Import", "DataFolder"), ("Import", "ProcessingFolder"), ("Import", "Extension"),
                          ("Import", "DosePerAngstromFrame"), ("Tomo", "DimensionsX"), ("Tomo", "DimensionsY"), ("Tomo", "DimensionsZ")):
        a, b = existing.get(section, name), wanted.get(section, name)
        if a != b:
            out.append(f"{section}/{name}: {a!r} vs {b!r}")
    return out


def ensure_settings(root: Path, wanted: WarpSettings) -> str:
    path = root / SETTINGS
    if path.exists():
        conflicts = _conflicts(WarpSettings.from_file(str(path)), wanted)
        if conflicts:
            raise ValueError(f"{path} disagrees with this series (one settings per Warp project):\n  " + "\n  ".join(conflicts))
        return "settings reused"
    root.mkdir(parents=True, exist_ok=True)
    wanted.to_file(str(path))
    return "settings written"


def _resolve_doc_path(p: Optional[str], doc_dir: Path) -> Optional[Path]:
    if not p:
        return None
    q = Path(p)
    return q if q.is_absolute() else (doc_dir / q)


def cets_to_warp(
    region,
    res: Resolver,
    sr: SeriesReport,
    *,
    root: Path,
    doc_dir: Path,
    companion: Optional[Companion],
    alignment_selector=None,
    tomogram_selector: Optional[str] = None,
    overwrite: bool = False,
) -> Dict[str, Path]:
    cets_alignment = select_alignment(region, alignment_selector)
    ts = find_by_id(region.tilt_series, cets_alignment.tilt_series_id, "tilt series")
    stem = ts.id
    aln_name = alignment_name_of(cets_alignment)
    comp_aln = companion.alignment(ts.id, aln_name) if companion else None
    comp_ts = companion.tilt_series.get(ts.id) if companion else None

    tomo = select_tomogram(region, cets_alignment, tomogram_selector, companion)
    reference = ReferenceVolume.from_tomogram(tomo)
    native = comp_aln.native_volume_dimension_a if comp_aln and comp_aln.native_volume_dimension_a else None
    hub: Alignment = alignment_from_cets(
        cets_alignment, tilt_series=ts, reference=reference, target_frame=FRAME_CONVENTIONS["WARP"],
        native_dimension_a=native, format_="WARP",
    )
    res.resolve("reference_tomogram", discovered=tomo.id, note="companion" if comp_aln and comp_aln.reference_tomogram_id else "region")

    images = sorted(ts.images or [], key=lambda im: im.section)
    n_rows = len(images)
    if [im.section for im in images] != list(range(n_rows)):
        raise ValueError(f"{stem}: tilt image sections are not 0..{n_rows - 1}")
    fr = image_frame(images[0])
    pix = fr.isotropic_spacing
    width, height = fr.size_px
    aligned = {p.z_index for p in hub.per_section_alignment_parameters}
    comp_images = comp_ts.images if comp_ts else {}

    acq_pix = res.require("pix", companion=(comp_ts.pixel_size_acquisition_a if comp_ts else None), discovered=pix,
                          note="acquisition pixel (settings Import/PixelSize); tilt-image pixel from the document")
    voltage = res.require("voltage", companion=(comp_ts.voltage_kv if comp_ts else None))
    cs = res.require("cs", companion=(comp_ts.cs_mm if comp_ts else None))
    amp = res.require("amp_contrast", companion=(comp_ts.amplitude_contrast if comp_ts else None))
    exposures = [comp_images[im.id].exposure_dose if im.id in comp_images else None for im in images]
    dose_per_tilt = res.require(
        "dose_per_tilt",
        companion=(comp_ts.dose_rate if comp_ts and comp_ts.dose_rate else (exposures[0] if all(e is not None for e in exposures) and len(set(exposures)) == 1 else None)),
        note="Warp holds one DosePerAngstromFrame per project",
    )
    angles_inverted = bool(res.optional("angles_inverted", companion=(comp_ts.are_angles_inverted if comp_ts else None), absent=False))

    # rows: nominal stage angle (Warp stores -stage), accumulated dose, movie path
    dark_angles = {im.section: float(im.nominal_tilt_angle or 0.0) for im in images if im.section not in aligned}
    doses = {im.section: float(im.accumulated_dose or 0.0) for im in images}
    if any(im.accumulated_dose is None for im in images):
        sr.warnings.append("some tilt images have no accumulated_dose: <Dose> written as 0 for them")
    stacks = {}
    if region.movie_stack_collection:
        for series in region.movie_stack_collection.movie_stacks or []:
            for st in series.stacks or []:
                stacks[st.id] = st
    frames_dir = Path(res.value("frames_dir", default=str(root / "frames")))
    tomostar_dir = root / TOMOSTAR_DIR
    movie_paths: Dict[int, str] = {}
    for im in images:
        name = None
        if im.movie_stack_id and im.movie_stack_id in stacks and stacks[im.movie_stack_id].path:
            name = Path(stacks[im.movie_stack_id].path.replace("\\", "/")).name
        elif im.id in comp_images and comp_images[im.id].frame_name:
            name = comp_images[im.id].frame_name
        if name is None:
            name = f"{stem}_{im.section:03d}.mrc"
        movie_paths[im.section] = os.path.relpath((frames_dir / name).resolve(), tomostar_dir.resolve())
    if any(im.movie_stack_id is None for im in images) and not comp_images:
        sr.warnings.append("no movie stacks / frame names in the document: MoviePath entries synthesised as <stem>_<section>.mrc")

    warp = hub.to_warp(pixel_size_a=pix, image_size_px=(width, height), n_rows=n_rows, dark_angles=dark_angles,
                       doses=doses, movie_paths=movie_paths)
    warp.are_angles_inverted = angles_inverted
    warp.data_directory = str(tomostar_dir.resolve())
    warp.ctf_params = {"PixelSize": f"{pix:.9g}", "Voltage": f"{voltage:g}", "Cs": f"{cs:g}", "Amplitude": f"{amp:g}"}
    ctfs = [im.ctf_metadata for im in images]
    if all(c is not None and c.defocus_u is not None for c in ctfs) and not res.optional("no_ctf", absent=False):
        for e, c in zip(warp.entries, ctfs, strict=True):
            e.defocus_um, e.defocus_delta_um, e.defocus_angle_deg, e.phase_shift_pi = cets_ctf.to_warp_values(c)
    elif any(c is not None for c in ctfs):
        sr.warnings.append("CTF metadata incomplete across tilt images: no CTF grids written")

    # operator gate: the written fields rebuild the CETS rotations
    check = Alignment.from_warp(warp)
    worst = 0.0
    for pa, p in zip(cets_alignment.projection_alignments, check.per_section_alignment_parameters, strict=True):
        r_cets, _ = fold_projection(pa)
        d = p.rotation_matrix() - r_cets
        worst = max(worst, float(np.degrees(np.arcsin(min(1.0, np.linalg.norm(d, 2) / 2.0)))))
    sr.gates.append(Gate("operators_rebuilt", worst <= OPERATOR_TOL_DEG, value=worst, expected=f"<= {OPERATOR_TOL_DEG} deg"))

    # settings (one per project)
    vol_a = np.array(reference.extent_a if native is None else [native["x"], native["y"], native["z"]])
    tomo_dims_px = [int(round(v / acq_pix)) for v in vol_a]
    settings = WarpSettings.create(
        pixel_size_a=float(acq_pix), exposure_per_tilt=float(dose_per_tilt), tomo_dims_px=tomo_dims_px,
        voltage_kv=float(voltage), cs_mm=float(cs), amplitude_contrast=float(amp),
    )
    note = ensure_settings(root, settings)
    sr.gates.append(Gate("settings", True, note=note))

    # tomostar
    rows = [
        TomostarRow(
            movie_name=movie_paths[im.section],
            angle_tilt=float(warp.entries[im.section].tilt_angle),
            axis_angle=float(hub.get_median_tilt_axis()),
            dose=doses[im.section],
        )
        for im in images
    ]
    tomostar = WarpTomostar(rows=rows)

    outputs = {"xml": root / PROCESSING_DIR / f"{stem}.xml", "tomostar": tomostar_dir / f"{stem}.tomostar", "settings": root / SETTINGS}
    for key in ("xml", "tomostar"):
        if outputs[key].exists() and not overwrite:
            raise FileExistsError(f"{outputs[key]} exists (use --overwrite)")
    tomostar_dir.mkdir(parents=True, exist_ok=True)
    (root / PROCESSING_DIR).mkdir(parents=True, exist_ok=True)
    tomostar.to_file(str(outputs["tomostar"]))
    warp.to_file(str(outputs["xml"]))
    reread = WarpAlignment.from_file(outputs["xml"], pixel_size_a=pix)
    sr.gates.append(Gate("rows_match_tomostar", reread.n_tilts == tomostar.n_rows == n_rows, value=(reread.n_tilts, tomostar.n_rows), expected=n_rows))

    missing = [movie_paths[im.section] for im in images if not (tomostar_dir / movie_paths[im.section]).exists()]
    if missing:
        sr.hints.append(f"place the movies (or their average/<name>.mrc) at the MoviePath targets, e.g. {tomostar_dir / missing[0]}")
    sr.hints.append(f"WarpTools ts_reconstruct --settings {root / SETTINGS} --angpix <voxel> --perdevice 1")
    sr.outputs.update({k: str(v) for k, v in outputs.items()})
    sr.provenance = res.provenance()
    return outputs
