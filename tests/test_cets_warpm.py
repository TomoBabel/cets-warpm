import json
import re
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner
from conftest import IMG, PIX, STEM, N, make_project
from cryoet_alignment.io.cets.companion import Companion
from cryoet_alignment.io.cets.entities import load_dataset
from cryoet_alignment.io.warp import WarpAlignment, WarpSettings, WarpTomostar

from cets_warpm.cli import main


def _run(args, expect_ok=True):
    result = CliRunner().invoke(main, args, catch_exceptions=False)
    if expect_ok:
        assert result.exit_code == 0, result.stdout + result.stderr
    return result


def test_to_cets_from_project_root(project, tmp_path):
    out = tmp_path / "o" / "warp.cets.json"
    r = _run(["to-cets", str(project), "-o", str(out)], expect_ok=False)
    assert r.exit_code == 1 and "spatially varying deformation grids" in r.stderr
    r = _run(["to-cets", str(project), "-o", str(out), "--drop-locals"])
    assert "pix = 1.54  [discovered]  (average/" in r.stdout
    assert "volume_dims_a = [6307.84, 6307.84, 3080.0]  [discovered]  (warp_tiltseries.settings" in r.stdout
    src = WarpAlignment.from_file(project / "warp_tiltseries" / f"{STEM}.xml", pixel_size_a=PIX, strict_dims=False)
    ds = load_dataset(out)
    region = ds.regions[0]
    ts = region.tilt_series[0]
    assert len(ts.images) == N and ts.images[0].width == IMG
    assert [im.nominal_tilt_angle for im in ts.images] == pytest.approx([-e.tilt_angle for e in src.entries])
    assert [im.accumulated_dose for im in ts.images] == pytest.approx([e.dose for e in src.entries])
    c = ts.images[0].ctf_metadata
    assert c.defocus_u == pytest.approx((src.entries[0].defocus_um + src.entries[0].defocus_delta_um / 2) * 1e4)
    assert c.defocus_angle == 17.0 and c.phase_shift == 0.0 and c.defocus_handedness is None
    assert ts.images[0].path.endswith(f"average/{Path(src.entries[0].movie_path).stem}.mrc")
    assert region.movie_stack_collection.movie_stacks[0].stacks[0].path.endswith(".eer")
    ids = [t.id for t in region.tomograms]
    assert ids == [f"{STEM}_volume", f"{STEM}_tomo_98.560"]
    assert (region.tomograms[0].width, region.tomograms[0].depth) == (4096, 2000)
    assert region.tomograms[1].coordinate_transformations[0].sequence[1].scale == pytest.approx([98.56] * 3)
    comp = Companion.load_for(out)
    tsc = comp.tilt_series[STEM]
    assert (tsc.voltage_kv, tsc.cs_mm, tsc.amplitude_contrast, tsc.dose_rate) == (300.0, 2.7, 0.07, 3.87)
    assert tsc.pixel_size_acquisition_a == 1.54 and tsc.pixel_size_ctf_a == 1.54 and tsc.are_angles_inverted is False
    assert tsc.tilt_axis_nominal_deg == pytest.approx(83.939545)
    order = [tsc.images[f"{STEM}_{i}"].acquisition_index_1b for i in range(N)]
    doses = [e.dose for e in src.entries]
    assert sorted(order) == list(range(1, N + 1)) and all(
        (doses[i] < doses[j]) == (order[i] < order[j]) for i in range(N) for j in range(N) if i != j
    )
    assert comp.alignments[0].reference_tomogram_id == f"{STEM}_volume" and comp.alignments[0].tomogram_ids == ids
    assert "GridMovementX" in comp.alignments[0].dropped[0]


def test_roundtrip_xml_fields_and_operators(project, tmp_path):
    out = tmp_path / "o" / "warp.cets.json"
    _run(["to-cets", str(project), "-o", str(out), "--drop-locals"])
    back = tmp_path / "back"
    r = _run(["from-cets", str(out), "-o", str(back)])
    assert "[ok ] operators_rebuilt" in r.stdout and "settings written" in r.stdout
    src = WarpAlignment.from_file(project / "warp_tiltseries" / f"{STEM}.xml", pixel_size_a=PIX, strict_dims=False)
    got = WarpAlignment.from_file(back / "warp_tiltseries" / f"{STEM}.xml", pixel_size_a=PIX)
    for f in (
        "tilt_angle",
        "tilt_axis_angle",
        "tilt_axis_offset_x",
        "tilt_axis_offset_y",
        "dose",
        "defocus_um",
        "defocus_delta_um",
        "defocus_angle_deg",
        "phase_shift_pi",
    ):
        a = np.array([getattr(e, f) for e in src.entries])
        b = np.array([getattr(e, f) for e in got.entries])
        assert np.abs(a - b).max() < 1e-6, f
    assert got.is_rigid and got.has_ctf and got.volume_dimensions_physical == pytest.approx([6307.84, 6307.84, 3080.0])
    assert got.ctf_params["Voltage"] == "300" and got.ctf_params["Cs"] == "2.7"
    settings = WarpSettings.from_file(str(back / "warp_tiltseries.settings"))
    assert (
        settings.tomo_dims_px == [4096, 4096, 2000]
        and settings.exposure_per_tilt == 3.87
        and settings.pixel_size_a == 1.54
    )
    tomostar = WarpTomostar.from_file(str(back / "tomostar" / f"{STEM}.tomostar"))
    assert tomostar.n_rows == N and [r.angle_tilt for r in tomostar.rows] == pytest.approx(
        [e.tilt_angle for e in src.entries]
    )
    assert tomostar.rows[0].movie_name == f"../frames/{Path(src.entries[0].movie_path).name}"
    # a second series into the same project must agree with the settings
    r = _run(["from-cets", str(out), "-o", str(back), "--overwrite", "--dose-per-tilt", "2.0"], expect_ok=False)
    assert "disagrees with this series" in r.stderr


def test_use_tilt_false_rows_and_constant_grid(tmp_path):
    def edit(xml: str) -> str:
        xml = re.sub(r"(<UseTilt>True\nTrue\n)True", r"\1False", xml, count=1)  # row 2 unused
        # replace the varying movement grids by a single-node constant one (7 Å) and a zero one
        xml = re.sub(
            r"<GridMovementX.*?</GridMovementX>",
            '<GridMovementX Width="1" Height="1" Depth="1" MarginX="0" MarginY="0" MarginZ="0"><Node X="0" Y="0" Z="0" Value="7" /></GridMovementX>',
            xml,
            flags=re.S,
        )
        xml = re.sub(
            r"<GridMovementY.*?</GridMovementY>",
            '<GridMovementY Width="1" Height="1" Depth="1" MarginX="0" MarginY="0" MarginZ="0"><Node X="0" Y="0" Z="0" Value="0" /></GridMovementY>',
            xml,
            flags=re.S,
        )
        return xml

    project = make_project(tmp_path, xml_edit=edit)
    out = tmp_path / "o" / "w.cets.json"
    r = _run(["to-cets", str(project), "-o", str(out)])  # rigid now: no --drop-locals needed
    assert "[ok ] constant_grids_folded" in r.stdout
    ds = load_dataset(out)
    pas = ds.regions[0].alignments[0].projection_alignments
    assert len(pas) == N - 1 and f"{STEM}_2" not in [pa.tilt_image_id for pa in pas]
    src = WarpAlignment.from_file(project / "warp_tiltseries" / f"{STEM}.xml", pixel_size_a=PIX, strict_dims=False)
    assert src.entries[0].movement_x == 7.0
    back = tmp_path / "back"
    _run(["from-cets", str(out), "-o", str(back)])
    got = WarpAlignment.from_file(back / "warp_tiltseries" / f"{STEM}.xml", pixel_size_a=PIX)
    assert [e.use_tilt for e in got.entries] == [e.use_tilt for e in src.entries]
    assert got.entries[2].tilt_angle == pytest.approx(src.entries[2].tilt_angle)  # dark row keeps its angle
    # the folded 7 Å landed in AxisOffsetX with a zero movement grid: same projection, different bookkeeping
    assert [e.tilt_axis_offset_x for e in got.entries if e.use_tilt] == pytest.approx(
        [e.tilt_axis_offset_x - 7.0 for e in src.entries if e.use_tilt]
    )
    assert all(e.movement_x == 0.0 for e in got.entries)


def test_zero_dims_need_settings_and_source_enumeration(tmp_path):
    project = make_project(tmp_path, with_recon=False)
    xml = project / "warp_tiltseries" / f"{STEM}.xml"
    r = _run(
        ["to-cets", str(xml), "-o", str(tmp_path / "a.cets.json"), "--drop-locals"]
    )  # settings found next to the processing folder
    assert r.exit_code == 0
    (project / "warp_tiltseries.settings").unlink()
    r = _run(["to-cets", str(xml), "-o", str(tmp_path / "b.cets.json"), "--drop-locals"], expect_ok=False)
    assert "volume_dims_a is required" in r.stderr  # image dims/pix still come from the averages
    r = _run(
        [
            "to-cets",
            str(xml),
            "-o",
            str(tmp_path / "c.cets.json"),
            "--drop-locals",
            "--volume-px",
            "4096x4096x2000",
            "--pix",
            "1.54",
        ]
    )
    assert "volume_dims_a = [6307.84, 6307.84, 3080.0]  [discovered]  (--volume-px x pixel)" in r.stdout
    # M enumeration: a .source listing this tomostar and a missing one, and a .population pointing at it
    src = project / "warp_tiltseries" / "proj.source"
    src.write_text(
        '<?xml version="1.0"?><DataSource><Param Name="PixelSize" Value="1.54"/><Param Name="DimensionsX" Value="4096"/>'
        '<Param Name="DimensionsY" Value="4096"/><Param Name="DimensionsZ" Value="2000"/><Param Name="DosePerAngstromFrame" Value="-3.87"/>'
        f'<Files><File Hash="a" Name="{STEM}.tomostar"/><File Hash="b" Name="Other.tomostar"/></Files></DataSource>',
    )
    (project / "proj.population").write_text(
        '<?xml version="1.0"?><Population><Sources><Source GUID="x" Path="warp_tiltseries/proj.source"/></Sources></Population>'
    )
    r = _run(["to-cets", str(project / "proj.population"), "-o", str(tmp_path / "d.cets.json"), "--drop-locals"])
    assert (
        "Other.tomostar: listed in the .source" in r.stderr
        and "volume_dims_a = [6307.84, 6307.84, 3080.0]  [discovered]  (.source" in r.stdout
    )
    assert len(load_dataset(tmp_path / "d.cets.json").regions) == 1


def test_from_cets_requires_inputs_without_companion(project, tmp_path):
    out = tmp_path / "o" / "warp.cets.json"
    _run(["to-cets", str(project), "-o", str(out), "--drop-locals"])
    Companion.path_for(out).unlink()
    r = _run(["from-cets", str(out), "-o", str(tmp_path / "b"), "--tomogram", f"{STEM}_volume"], expect_ok=False)
    assert "voltage is required" in r.stderr
    cfg = tmp_path / "cets.yaml"
    cfg.write_text(
        "cets:\n  voltage: 300\n  cs: 2.7\n  amp_contrast: 0.07\n  dose_per_tilt: 3.87\n  frames_dir: movies\n"
    )
    r = _run(["from-cets", str(out), "-o", str(tmp_path / "b"), "--tomogram", f"{STEM}_volume", "--config", str(cfg)])
    assert "voltage = 300  [config]" in r.stdout and "defaulted" not in r.stderr
    got = WarpAlignment.from_file(tmp_path / "b" / "warp_tiltseries" / f"{STEM}.xml", pixel_size_a=PIX)
    assert got.entries[0].movie_path.startswith("../movies/") or "movies" in got.entries[0].movie_path
    assert got.is_rigid
    report = json.loads((tmp_path / "b" / "cets_warpm.report.json").read_text())
    assert report["ok"]
