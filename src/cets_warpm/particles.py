"""Particle stars ⇄ CETS annotations for Warp / M projects.

Import (``to-cets --particles``): every star (Warp import, M species import or RELION 5 flavour, see
``cryoet_alignment.io.cets.particles_star``) is read into corner-anchored Å positions and particle→tomogram
matrices, its rows are grouped by tilt-series stem and bound to that region's reference tomogram
(``<stem>_volume``, Warp's own box, or ``--particles-tomogram``), and become one ``PointSet3D`` /
``PointMatrixSet3D`` per (series, star) in the tomogram's centred physical frame.

Export (``from-cets --particles-out``): point annotations bound to the exported alignment's reference tomogram
are written as one star per document (or per series) in the requested flavour, with the coordinates in pixels
of ``--coords-angpix`` (default: the tomogram's voxel size, warned) measured from the corner of the Warp box.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
from cryoet_alignment.io.cets.annotations import (
    POINT_KINDS,
    ResolvedPoints,
    annotation_id,
    annotation_kind,
    annotation_points,
    point_set_entity,
    select_annotations,
    tomogram_frame,
)
from cryoet_alignment.io.cets.cli_support import Gate, SeriesReport
from cryoet_alignment.io.cets.companion import AnnotationCompanion, Companion
from cryoet_alignment.io.cets.config import Resolver
from cryoet_alignment.io.cets.frames import find_by_id
from cryoet_alignment.io.cets.particles_star import (
    FLAVOURS,
    ParticleTable,
    StarRows,
    read_particle_star,
    safe_stem,
    validate_series_names,
    write_particle_star,
)

SOURCE_TOOL = {"warp": "Warp", "m": "M", "relion5": "RELION 5"}


# --------------------------------------------------------------------------- import


@dataclass
class StarImport:
    """Annotations + companions produced from one star, keyed by region (tilt-series stem)."""

    path: Path
    table: ParticleTable
    annotations: Dict[str, list] = field(default_factory=dict)
    companions: Dict[str, AnnotationCompanion] = field(default_factory=dict)


def _reference_tomogram(region, selector: Optional[str]):
    if selector:
        return find_by_id(region.tomograms, selector, "tomogram")
    want = f"{region.id}_volume"
    for t in region.tomograms or []:
        if t.id == want:
            return t
    if len(region.tomograms or []) == 1:
        return region.tomograms[0]
    raise ValueError(
        f"region {region.id!r}: no reference tomogram {want!r} and {len(region.tomograms or [])} tomograms; "
        "pass --particles-tomogram",
    )


def import_star(
    path: Path,
    regions: Dict[str, Any],
    res: Resolver,
    sr: SeriesReport,
    *,
    tomogram_selector: Optional[str] = None,
) -> StarImport:
    """Read one star and bind its rows to the regions of the document (``regions`` keyed by tilt-series stem)."""
    flavour = res.optional("star_flavour", absent="auto")
    if flavour not in ("auto", *FLAVOURS):
        raise ValueError(f"--star-flavour {flavour!r}: expected auto|{'|'.join(FLAVOURS)}")
    coords_angpix = res.optional("coords_angpix")
    angpix_shifts = res.optional("angpix_shifts")
    table = read_particle_star(path, flavour, coords_angpix=coords_angpix, angpix_shifts=angpix_shifts)
    res.resolve(
        "coords_angpix_used",
        discovered=table.coords_angpix,
        note=f"{table.coords_angpix_source} ({table.flavour} flavour)",
    )
    unknown = validate_series_names(table, list(regions))
    skip_unknown = bool(res.optional("skip_unknown_series", absent=False))
    if unknown:
        if not skip_unknown:
            raise ValueError(
                f"{path.name}: {len(unknown)} series name(s) are not in the document: {unknown[:5]}"
                f"{' ...' if len(unknown) > 5 else ''} (pass --skip-unknown-series to drop those rows)",
            )
        n_drop = int(sum(1 for s in table.series if s not in regions))
        sr.warnings.append(f"{path.name}: {n_drop} row(s) of unknown series dropped: {unknown[:5]}")
    sr.gates.append(
        Gate(
            "rows_bound",
            not unknown or skip_unknown,
            value=table.n - int(sum(1 for s in table.series if s not in regions)),
            expected=table.n,
        )
    )

    imp = StarImport(path=path, table=table)
    key = safe_stem(path.name)
    n_inside = 0
    n_total = 0
    for stem in sorted(set(table.series)):
        if stem not in regions:
            continue
        region = regions[stem]
        rows = table.rows_for(stem)
        tomo = _reference_tomogram(region, tomogram_selector)
        frame = tomogram_frame(tomo)
        corner = table.positions_corner_a[rows]
        if table.centred_input:  # RELION 5 centred coordinates are relative to the float centre of the box
            corner = corner + frame.extent_a / 2.0
        inside = frame.inside(corner)
        n_inside += int(inside.sum())
        n_total += len(rows)
        pts = frame.corner_to_cets(corner)
        mats = table.matrices[rows] if table.matrices is not None else None
        ann_id = annotation_id(tomo.id, key)
        existing = {a.id for a in region.annotations or []}
        if ann_id in existing:
            raise ValueError(f"region {region.id!r}: annotation {ann_id!r} already exists (star stems must be unique)")
        ann = point_set_entity(annotation_id=ann_id, tomogram_id=tomo.id, points_a=pts, matrices=mats, name=path.stem)
        imp.annotations.setdefault(stem, []).append(ann)
        imp.companions[ann_id] = AnnotationCompanion(
            kind="oriented_points" if mats is not None else "points",
            tomogram_id=tomo.id,
            source_tool=SOURCE_TOOL[table.flavour],
            source_ref=str(path),
            name=path.stem,
            flavour=table.flavour,
            coords_angpix_a=table.coords_angpix,
            angpix_shifts_a=table.angpix_shifts,
            series_name=table.series_raw[int(rows[0])],
            columns={c: [v[i] for i in rows] for c, v in table.extra.items()},
            dropped=list(table.dropped),
        )
    if n_inside != n_total:
        sr.warnings.append(f"{path.name}: {n_total - n_inside} of {n_total} points lie outside the bound tomogram box")
    res.resolve(
        "points_inside_box", discovered=f"{n_inside}/{n_total}", note="corner-anchored box of the bound tomogram"
    )
    if table.matrices is None:
        sr.warnings.append(f"{path.name}: no rlnAngleRot/Tilt/Psi columns: written as PointSet3D (no orientations)")
    if table.dropped:
        sr.dropped.append(f"star columns not preserved: {table.dropped}")
    sr.provenance = res.provenance()
    return imp


# --------------------------------------------------------------------------- export


@dataclass
class ExportSelection:
    region_id: str
    stem: str
    resolved: ResolvedPoints
    companion: Optional[AnnotationCompanion]


def collect_for_export(
    region,
    reference_tomogram,
    companion: Optional[Companion],
    sr: SeriesReport,
    *,
    annotation_ids: Sequence[str] = (),
) -> List[ExportSelection]:
    """Point annotations of a region to export: those bound to the reference tomogram, or the ids requested."""
    if annotation_ids:
        anns = select_annotations(region, ids=list(annotation_ids))
    else:
        anns = select_annotations(region, tomogram_id=str(reference_tomogram.id), kinds=list(POINT_KINDS))
    out: List[ExportSelection] = []
    for ann in anns:
        kind = annotation_kind(ann)
        if kind not in POINT_KINDS:
            sr.warnings.append(f"annotation {ann.id!r} ({kind}) has no Warp/M/RELION particle representation - skipped")
            continue
        r = annotation_points(ann, region)
        if str(r.tomogram.id) != str(reference_tomogram.id):
            ref = tomogram_frame(reference_tomogram)
            if not np.allclose(r.frame.extent_a, ref.extent_a, atol=max(r.frame.spacing_a, ref.spacing_a)):
                raise ValueError(
                    f"annotation {ann.id!r} is bound to tomogram {r.tomogram.id!r} (extent {r.frame.extent_a.round(2).tolist()} Å), "
                    f"not to the exported reference {reference_tomogram.id!r} (extent {ref.extent_a.round(2).tolist()} Å)",
                )
            sr.warnings.append(
                f"annotation {ann.id!r} is bound to {r.tomogram.id!r}; converted through the shared centred frame of "
                f"{reference_tomogram.id!r} (same extent)",
            )
            # same physical extent: the centred frames coincide up to the floor(N/2) difference of the two grids
            shift = r.frame.corner_a - ref.corner_a  # move the CETS origin of the source grid onto the reference grid
            r = ResolvedPoints(r.annotation_id, reference_tomogram, ref, r.points_a + shift, r.matrices)
        out.append(
            ExportSelection(
                region_id=str(region.id),
                stem=str(region.id),
                resolved=r,
                companion=companion.annotations.get(str(ann.id)) if companion else None,
            ),
        )
    return out


def write_stars(
    selections: List[ExportSelection],
    res: Resolver,
    sr: SeriesReport,
    *,
    out_dir: Path,
    doc_stem: str,
    per_series: bool,
    settings_path: Optional[Path],
    overwrite: bool,
) -> Dict[str, Path]:
    """Write the selected annotations as particle star(s) and re-read them as a round-trip gate."""
    flavour = res.require("star_flavour")
    if flavour not in FLAVOURS:
        raise ValueError(f"--star-flavour {flavour!r}: expected {'|'.join(FLAVOURS)}")
    style = res.value("series_name_style", default="tomostar")
    if not selections:
        sr.warnings.append("no point annotations bound to the exported tomograms: no star written")
        return {}
    voxels = sorted({round(s.resolved.frame.spacing_a, 6) for s in selections})
    a = res.value(
        "coords_angpix",
        default=voxels[0],
        note="the bound tomogram's voxel size" + (f"; several: {voxels}" if len(voxels) > 1 else ""),
    )
    if len(voxels) > 1 and "coords_angpix" not in res.cli and not per_series:
        raise ValueError(
            f"tomogram voxel sizes differ between regions ({voxels}); pass --coords-angpix or --per-series"
        )
    extents = sorted({tuple(np.round(s.resolved.frame.extent_a, 3)) for s in selections})
    if flavour == "relion5" and len(extents) > 1 and not per_series:
        raise ValueError(f"box extents differ between regions ({extents}); relion5 needs --per-series")

    groups: Dict[str, List[ExportSelection]] = {}
    for s in selections:
        groups.setdefault(s.stem if per_series else doc_stem, []).append(s)
    written: Dict[str, Path] = {}
    n_points = 0
    worst_pos = 0.0
    worst_rot = 0.0
    for gname, sel in groups.items():
        series: List[str] = []
        pos: List[np.ndarray] = []
        mats: List[np.ndarray] = []
        extra: Dict[str, list] = {}
        any_mats = any(s.resolved.matrices is not None for s in sel)
        if any_mats and not all(s.resolved.matrices is not None for s in sel):
            sr.warnings.append(f"{gname}: mixing oriented and plain point sets; Eulers written as 0 for the plain ones")
        cols: Optional[set] = None
        for s in sel:
            n = s.resolved.n
            series += [s.stem] * n
            pos.append(s.resolved.points_corner_a)
            mats.append(
                s.resolved.matrices if s.resolved.matrices is not None else np.repeat(np.eye(3)[None], n, axis=0)
            )
            c: Dict[str, list] = s.companion.columns if s.companion else {}
            cols = set(c) if cols is None else (cols & set(c))
            for k, v in c.items():
                extra.setdefault(k, []).extend(v)
        # keep only columns every selection carries with the right length
        extra = {k: v for k, v in extra.items() if cols and k in cols and len(v) == len(series)}
        if not any_mats:
            sr.warnings.append(f"{gname}: point sets without orientations; rlnAngleRot/Tilt/Psi written as 0")
        rows = StarRows(
            series=series,
            positions_corner_a=np.concatenate(pos),
            matrices=np.concatenate(mats),
            extra=extra,
        )
        extent = sel[0].resolved.frame.extent_a
        path = out_dir / f"{gname}_{flavour}.star"
        if path.exists() and not overwrite:
            raise FileExistsError(f"{path} exists (use --overwrite)")
        write_particle_star(
            path,
            flavour,
            rows,
            coords_angpix=float(a),
            series_name_style=style,
            extent_a=extent,
            voltage_kv=res.optional("voltage"),
            cs_mm=res.optional("cs"),
            amplitude_contrast=res.optional("amp_contrast"),
        )
        written[gname] = path
        # round trip
        back = read_particle_star(path, flavour, coords_angpix=float(a) if flavour != "relion5" else None)
        got = back.positions_corner_a + (extent / 2.0 if back.centred_input else 0.0)
        worst_pos = max(worst_pos, float(np.abs(got - rows.positions_corner_a).max()))
        if back.matrices is not None:
            rel = np.einsum("nij,nkj->nik", back.matrices, rows.matrices)
            worst_rot = max(worst_rot, float(np.abs(rel - np.eye(3)).max()))
        n_points += len(series)
    sr.gates.append(
        Gate("star_roundtrip_positions", worst_pos < 1e-3, value=worst_pos, expected="< 1e-3 A (file precision)")
    )
    sr.gates.append(
        Gate("star_roundtrip_rotations", worst_rot < 1e-6, value=worst_rot, expected="< 1e-6 (|R_out R_in^T - I|)")
    )
    sr.outputs.update({f"star:{k}": str(v) for k, v in written.items()})
    first = next(iter(written.values()))
    settings = settings_path if settings_path else Path("<root>/warp_tiltseries.settings")
    if flavour == "warp":
        sr.hints.append(
            f"WarpTools ts_export_particles --settings {settings} --input_star {first} --coords_angpix {float(a):g} "
            f"--output_star {first.with_name(first.stem + '_export.star')} --box <px> --diameter <A> --2d",
        )
    elif flavour == "m":
        sr.hints.append(
            f"MTools create_species --population <project>.population --particles_relion {first} "
            f"--angpix_coords {float(a):g} --name <species> --diameter <A> --sym C1 --half1 <mrc> --half2 <mrc>",
        )
    else:
        sr.hints.append(
            f"RELION 5: particles star {first} (rlnCenteredCoordinate*Angst, rlnTomoName = <stem>{'.tomostar' if style == 'tomostar' else ''})"
        )
    sr.provenance = res.provenance()
    return written


__all__ = [
    "StarImport",
    "import_star",
    "ExportSelection",
    "collect_for_export",
    "write_stars",
    "SOURCE_TOOL",
]
