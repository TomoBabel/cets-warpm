"""A Warp 2.0 project assembled from real text fixtures (24jul16a Position_16_3: tilt-series XML with 4x4x31
movement grids and 1x1x31 CTF grids, settings, tomostar) plus synthetic frame-average / reconstruction MRC
headers (128 px images at 1.54 Å; a 64-voxel reconstruction)."""

import re
import shutil
from pathlib import Path

import mrcfile
import numpy as np
import pytest

DATA = Path(__file__).parent / "data"
STEM = "Position_16_3"
PIX = 1.54
IMG = 128
N = 31


def write_mrc(path: Path, shape_zyx, voxel: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with mrcfile.new(str(path), overwrite=True) as m:
        m.set_data(np.zeros(shape_zyx, dtype=np.float32))
        m.voxel_size = voxel


def make_project(root: Path, *, xml_edit=None, with_recon=True, with_averages=True) -> Path:
    root = root / "warp"
    (root / "tomostar").mkdir(parents=True, exist_ok=True)
    (root / "warp_tiltseries").mkdir(parents=True, exist_ok=True)
    shutil.copy(DATA / "warp_tiltseries.settings", root / "warp_tiltseries.settings")
    shutil.copy(DATA / f"{STEM}.tomostar", root / "tomostar" / f"{STEM}.tomostar")
    xml = (DATA / f"{STEM}.xml").read_text()
    if xml_edit is not None:
        xml = xml_edit(xml)
    (root / "warp_tiltseries" / f"{STEM}.xml").write_text(xml)
    if with_averages:
        for line in re.search(r"<MoviePath>(.*?)</MoviePath>", xml, re.S).group(1).strip().split("\n"):
            movie = (root / "tomostar" / line.strip()).resolve()
            write_mrc(movie.parent / "average" / (movie.stem + ".mrc"), (1, IMG, IMG), PIX)
    if with_recon:
        write_mrc(root / "warp_tiltseries" / "reconstruction" / f"{STEM}_98.56Apx.mrc", (32, 64, 64), 98.56)
    return root


@pytest.fixture
def project(tmp_path) -> Path:
    return make_project(tmp_path)
