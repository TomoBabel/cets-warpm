"""Particle stars <-> CETS annotations through the cets-warpm CLI: the three flavours round-trip, the M import
formula against M's own species table, unknown series, config provenance, export selection rules."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import starfile
from click.testing import CliRunner
from conftest import PIX, STEM, make_project
from cryoet_alignment.io.cets.annotations import annotation_points, tomogram_frame
from cryoet_alignment.io.cets.companion import Companion
from cryoet_alignment.io.cets.entities import load_dataset
from cryoet_alignment.io.cets.euler import zyz_to_matrices
from cryoet_alignment.io.cets.particles_star import read_particle_star

from cets_warpm.cli import main

DATA = Path(__file__).parent / "data" / "particles"


def _run(args, expect_ok=True):
    result = CliRunner().invoke(main, args, catch_exceptions=False)
    if expect_ok:
        assert result.exit_code == 0, result.stdout + result.stderr
    return result


def _write_picks(path: Path, name=f"{STEM}.tomostar", eulers=True, extra=True):
    df = pd.DataFrame(
        {
            "rlnMicrographName": [name] * 3,
            "rlnCoordinateX": [1000.5, 1500.0, 100.0],
            "rlnCoordinateY": [2000.25, 1000.0, 400.0],
            "rlnCoordinateZ": [300.0, 700.0, 150.0],
        }
    )
    if eulers:
        df["rlnAngleRot"], df["rlnAngleTilt"], df["rlnAnglePsi"] = [10.0, -50.0, 0.0], [20.0, 100.0, 90.0], [30.0, 120.0, 0.0]
    if extra:
        df["rlnRandomSubset"] = [1, 2, 1]
        df["rlnSomethingElse"] = [0.1, 0.2, 0.3]
    starfile.write({"particles": df}, path)
    return df


@pytest.mark.parametrize("flavour", ["warp", "m", "relion5"])
def test_star_import_export_roundtrip(tmp_path, flavour):
    project = make_project(tmp_path)
    picks = tmp_path / "picks.star"
    df = _write_picks(picks)
    out = tmp_path / "cets" / "warp.cets.json"
    r = _run(["to-cets", str(project), "-o", str(out), "--drop-locals", "--particles", str(picks), "--coords-angpix", "4.0"])
    assert "coords_angpix_used = 4.0  [discovered]  (explicit (warp flavour))" in r.stdout
    assert "[ok ] rows_bound value=3 expected=3" in r.stdout
    assert "points lie outside" not in r.stdout + r.stderr or "of 3 points" in r.stderr  # depends on the box
    ds = load_dataset(out)
    region = ds.regions[0]
    assert len(region.annotations) == 1 and region.annotations[0].id == f"{STEM}_volume_ann_picks"
    comp = Companion.load_for(out)
    ac = comp.annotations[region.annotations[0].id]
    assert ac.flavour == "warp" and ac.coords_angpix_a == 4.0 and ac.columns["rlnRandomSubset"] == [1, 2, 1]
    assert "rlnSomethingElse" in ac.dropped
    # the CETS points are the corner-anchored positions minus the CETS centre of the Warp box
    resolved = annotation_points(region.annotations[0], region)
    frame = tomogram_frame(resolved.tomogram)
    expect_corner = df[["rlnCoordinateX", "rlnCoordinateY", "rlnCoordinateZ"]].to_numpy(float) * 4.0
    assert np.allclose(resolved.points_corner_a, expect_corner, atol=1e-6)
    assert np.allclose(resolved.matrices, zyz_to_matrices(df[["rlnAngleRot", "rlnAngleTilt", "rlnAnglePsi"]].to_numpy(float)), atol=1e-12)
    assert np.allclose(frame.corner_a, [float(n // 2) * PIX for n in frame.size_px])

    root = tmp_path / "warp_out"
    r = _run(
        [
            "from-cets", str(out), "-o", str(root), "--particles-out", str(root / "particles"), "--star-flavour", flavour,
            "--coords-angpix", "2.5", "--frames-dir", str(root / "frames"), "--voltage", "300", "--cs", "2.7", "--amp-contrast", "0.07",
        ]
    )
    assert "[ok ] star_roundtrip_positions" in r.stdout and "[ok ] star_roundtrip_rotations" in r.stdout
    star = root / "particles" / f"warp_{flavour}.star"
    assert star.exists()
    back = read_particle_star(star, flavour, coords_angpix=2.5 if flavour != "relion5" else None)
    corner = back.positions_corner_a + (frame.extent_a / 2.0 if back.centred_input else 0.0)
    assert np.allclose(corner, expect_corner, atol=1e-4)
    assert np.abs(np.einsum("nij,nkj->nik", back.matrices, resolved.matrices) - np.eye(3)).max() < 1e-7
    assert back.extra["rlnRandomSubset"] == [1, 2, 1]
    raw = starfile.read(star, always_dict=True)
    if flavour == "relion5":
        assert raw["particles"]["rlnTomoName"].iloc[0] == f"{STEM}.tomostar"
        assert raw["optics"]["rlnTomoTiltSeriesPixelSize"].iloc[0] == pytest.approx(2.5)
        assert "RELION 5" in r.stdout
    else:
        assert raw["particles"]["rlnMicrographName"].iloc[0] == f"{STEM}.tomostar"
        assert raw["particles"]["rlnCoordinateX"].iloc[0] == pytest.approx(1000.5 * 4.0 / 2.5, abs=1e-5)
        hint = "ts_export_particles" if flavour == "warp" else "create_species"
        assert hint in r.stdout and "--coords_angpix 2.5" in r.stdout or "--angpix_coords 2.5" in r.stdout


def test_m_import_star_binds_with_the_m_pixel_chain(tmp_path):
    """A RELION-refined star read as the M species import: rlnImagePixelSize (4.0), not the tilt-series pixel."""
    project = make_project(tmp_path)
    star = DATA / "m_import_run_data.star"
    out = tmp_path / "cets" / "warp.cets.json"
    r = _run(["to-cets", str(project), "-o", str(out), "--drop-locals", "--particles", str(star), "--star-flavour", "m"], expect_ok=False)
    assert r.exit_code == 1 and "series name(s) are not in the document" in r.stderr
    r = _run(["to-cets", str(project), "-o", str(out), "--drop-locals", "--particles", str(star), "--star-flavour", "m", "--skip-unknown-series", "--overwrite"])
    assert "coords_angpix_used = 4.0  [discovered]  (rlnImagePixelSize (m flavour))" in r.stdout
    assert "row(s) of unknown series dropped" in r.stderr
    assert load_dataset(out).regions[0].annotations == []  # none of the excerpt's series is in this project
    # the reader itself against M's table (the CLI just routes it)
    t = read_particle_star(star, "m")
    mt = starfile.read(DATA / "m_species_particles.star")
    assert np.abs(mt[["wrpCoordinateX1", "wrpCoordinateY1", "wrpCoordinateZ1"]].to_numpy(float) - t.positions_corner_a).max() < 1e-3


def test_star_options_from_config_and_no_orientations(tmp_path):
    project = make_project(tmp_path)
    picks = tmp_path / "picks.star"
    _write_picks(picks, eulers=False, extra=False)
    cfg = tmp_path / "cets.yaml"
    cfg.write_text("cets-warpm:\n  to-cets: {coords_angpix: 3.0, star_flavour: warp}\n  from-cets: {coords_angpix: 3.0, series_name_style: stem}\n")
    out = tmp_path / "cets" / "warp.cets.json"
    r = _run(["to-cets", str(project), "-o", str(out), "--drop-locals", "--particles", str(picks), "--config", str(cfg)])
    assert "coords_angpix = 3.0  [config]" in r.stdout and "star_flavour = 'warp'  [config]" in r.stdout
    assert "written as PointSet3D (no orientations)" in r.stderr
    root = tmp_path / "warp_out"
    r = _run(["from-cets", str(out), "-o", str(root), "--particles-out", str(root / "p"), "--star-flavour", "warp", "--config", str(cfg), "--frames-dir", str(root / "frames")])
    assert "coords_angpix = 3.0  [config]" in r.stdout and "WARNING: coords_angpix defaulted" not in r.stderr
    assert "rlnAngleRot/Tilt/Psi written as 0" in r.stderr
    raw = starfile.read(root / "p" / "warp_warp.star")
    assert raw["rlnMicrographName"].iloc[0] == STEM and set(raw["rlnAngleRot"]) == {0.0}


def test_from_cets_without_particles_out_writes_no_star_and_flavour_required(tmp_path):
    project = make_project(tmp_path)
    picks = tmp_path / "picks.star"
    _write_picks(picks)
    out = tmp_path / "cets" / "warp.cets.json"
    _run(["to-cets", str(project), "-o", str(out), "--drop-locals", "--particles", str(picks), "--coords-angpix", "4.0"])
    root = tmp_path / "warp_out"
    r = _run(["from-cets", str(out), "-o", str(root), "--frames-dir", str(root / "frames")])
    assert not (root / "particles").exists() and "star_roundtrip" not in r.stdout
    r = _run(["from-cets", str(out), "-o", str(root), "--particles-out", str(root / "p"), "--frames-dir", str(root / "frames"), "--overwrite"], expect_ok=False)
    assert r.exit_code == 1 and "star_flavour is required" in r.stderr
