"""X1: a real AreTomo3 .aln encoded into CETS by the codec, written as a Warp project by cets-warpm, loaded by
arewarpion's INDEPENDENT Warp model, must project 3D points like arewarpion's AreTomo3 model of the source
.aln (float32 warpylib: 1e-2 Å)."""

from pathlib import Path

import numpy as np
import pytest

torch = pytest.importorskip("torch")
warp_xml = pytest.importorskip("arewarpion.io.warp_xml")
warp_ts = pytest.importorskip("arewarpion.models.warp_ts")
aretomo_ts = pytest.importorskip("arewarpion.models.aretomo_ts")

from cryoet_alignment.io.aretomo3 import AreTomo3ALN  # noqa: E402
from cryoet_alignment.io.cets.alignment import ReferenceVolume, alignment_to_cets  # noqa: E402
from cryoet_alignment.io.cets.cli_support import SeriesReport  # noqa: E402
from cryoet_alignment.io.cets.config import Resolver  # noqa: E402
from cryoet_alignment.io.cets.entities import region_entity, tilt_series_entity, tomogram_entity  # noqa: E402
from cryoet_alignment.io.cets.frames import FRAME_CONVENTIONS, image_frame  # noqa: E402
from cryoet_alignment.io.cryoet_data_portal import Alignment  # noqa: E402

from cets_warpm.from_cets import cets_to_warp  # noqa: E402

ALN = Path("/hpc/projects/group.czii/utz.ermel/repos/arewarpo/testdata_runs/at3_24jul16a/outB/24jul16a_Position_16_3.aln")


@pytest.mark.skipif(not ALN.exists(), reason="real AreTomo3 run not available")
def test_x1_aretomo3_to_cets_to_warp_projects_identically(tmp_path):
    aln = AreTomo3ALN.from_file(ALN)
    s, vol = 1.54, (4096, 4096, 2000)
    hub = Alignment.from_aretomo3(aln, vol_size_px=vol, pixel_size_a=s)
    ts = tilt_series_entity(tilt_series_id="TS", path=None, width=4096, height=4096, pixel_size_a=s,
                            nominal_angles=[t - aln.AlphaOffset for t in aln.raw_tilts()], doses=[3.0 * i for i in range(aln.n_raw)])
    ref = tomogram_entity(tomogram_id="TS_volume", path=None, size_px=vol, voxel_size_a=s, tilt_series_id="TS")
    cets = alignment_to_cets(hub, tilt_series_id="TS", alignment_name="aretomo3", image=image_frame(ts.images[0]),
                             reference=ReferenceVolume.from_tomogram(ref), frame=FRAME_CONVENTIONS["ARETOMO3"])
    region = region_entity(region_id="TS", tilt_series=[ts], alignments=[cets], tomograms=[ref])
    sr = SeriesReport("TS")
    res = Resolver("cets-warpm", "from-cets", cli={"voltage": 300.0, "cs": 2.7, "amp_contrast": 0.07, "dose_per_tilt": 3.0, "pix": s}, warn=sr.warnings.append)
    outputs = cets_to_warp(region, res, sr, root=tmp_path / "warp", doc_dir=tmp_path, companion=None)
    assert all(g.passed for g in sr.gates), [g for g in sr.gates if not g.passed]

    a3 = aretomo_ts.AretomoTsModel(
        rot_deg=torch.tensor([g.rot for g in aln.GlobalAlignments], dtype=torch.float64),
        tilt_deg=torch.tensor([g.tilt for g in aln.GlobalAlignments], dtype=torch.float64),
        shifts_px=torch.tensor([[g.tx, g.ty] for g in aln.GlobalAlignments], dtype=torch.float64),
        raw_size_px=(4096, 4096), pixel_size_a=s, volume_dims_a=tuple(v * s for v in vol),
    )
    wm = warp_ts.WarpTiltSeriesModel(warp_xml.load_warp_tiltseries(outputs["xml"]).ts)
    rng = np.random.default_rng(1)
    pts = torch.tensor(rng.uniform(0.15, 0.85, size=(300, 3)) * np.array(vol) * s, dtype=torch.float64)
    xa, _ = a3.project_volume_global(pts)
    xw, _ = wm.project_volume_global(pts)
    d = (xa - xw).abs().max().item()
    assert d < 1e-2, f"AreTomo3 model vs Warp model of the CETS-written XML: {d:.4f} Å"
