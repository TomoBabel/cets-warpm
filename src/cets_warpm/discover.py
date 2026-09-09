"""Locate the files of a Warp 2.0 tilt-series project (or an M population/source pointing at one) and read
what each tilt series needs, reporting every value with its source.

Pixel sizes are distinct quantities and are kept apart: ``Import/PixelSize`` (acquisition),
``PixelSize × 2^BinTimes`` (stored frame averages = tilt images; verified against the average header),
``<CTF PixelSize>`` (the sampling Warp used for CTF fitting), and the reconstruction voxel (from the
reconstruction file). Volume box: XML ``VolumeDimensionsAngstrom`` when non-zero, else
``Tomo/Dimensions × Import/PixelSize`` from the settings / .source (re-saved XMLs carry "0, 0, 0").
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from xml.etree import ElementTree

import mrcfile
from cryoet_alignment.io.warp import WarpAlignment, WarpSettings, WarpTomostar

SETTINGS_NAME = "warp_tiltseries.settings"


@dataclass
class Found:
    value: Any
    source: str


@dataclass
class WarpSeries:
    stem: str
    xml_path: Path
    settings_path: Optional[Path] = None
    tomostar_path: Optional[Path] = None
    tomostar_dir: Optional[Path] = None
    reconstructions: List[Path] = field(default_factory=list)
    found: Dict[str, Found] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    settings: Optional[WarpSettings] = None
    tomostar: Optional[WarpTomostar] = None

    def get(self, key: str, default=None):
        f = self.found.get(key)
        return default if f is None else f.value

    def source(self, key: str) -> str:
        f = self.found.get(key)
        return "" if f is None else f.source

    def add(self, key: str, value: Any, source: str, override: bool = False) -> None:
        if value is None:
            return
        if key not in self.found or override:
            self.found[key] = Found(value, source)


def mrc_header(path: Path) -> dict:
    with mrcfile.open(str(path), header_only=True, permissive=True) as m:
        h = m.header
        vs = m.voxel_size
        return {"nx": int(h.nx), "ny": int(h.ny), "nz": int(h.nz), "voxel": (float(vs.x), float(vs.y), float(vs.z))}


def _read_source(path: Path) -> Tuple[Dict[str, str], List[str]]:
    """M ``.source``: ``{Param name: value}`` and the tomostar file names (relative to the .source folder)."""
    root = ElementTree.parse(path).getroot()
    params = {p.get("Name"): p.get("Value", "") for p in root.findall("Param")}
    files = [f.get("Name") for f in root.findall("Files/File")]
    return params, files


def _read_population(path: Path) -> List[Path]:
    root = ElementTree.parse(path).getroot()
    return [(path.parent / s.get("Path")).resolve() for s in root.findall("Sources/Source")]


def enumerate_xmls(tokens: List[str], settings_override: Optional[Path] = None) -> List[Tuple[Path, Optional[Path], Optional[dict]]]:
    """Expand CLI sources to ``(xml_path, settings_path, source_params)`` triples.

    Tokens: a project root (``<root>/warp_tiltseries.settings``), a ``.settings`` file, ``.xml`` files or
    directories of them, an M ``.source`` (its tomostars → ``<folder>/<stem>.xml``) or ``.population``.
    """
    out = []
    seen = set()

    def push(xml: Path, settings: Optional[Path], src: Optional[dict]):
        key = xml.resolve()
        if key not in seen:
            seen.add(key)
            out.append((xml, settings, src))

    def from_settings(settings: Path):
        st = WarpSettings.from_file(str(settings))
        proc = settings.parent / (st.processing_folder or "warp_tiltseries")
        xmls = sorted(proc.glob("*.xml"))
        if not xmls:
            raise FileNotFoundError(f"{proc}: no tilt-series XML files")
        for x in xmls:
            push(x, settings, None)

    def from_source(source: Path):
        params, files = _read_source(source)
        settings = settings_override or _find_settings(source.parent)
        for name in files:
            xml = source.parent / (Path(name).stem + ".xml")
            if xml.exists():
                push(xml, settings, params)
            else:
                out.append((xml, settings, {"__missing__": name}))

    for tok in tokens:
        p = Path(tok)
        if p.is_dir():
            if (p / SETTINGS_NAME).exists():
                from_settings(p / SETTINGS_NAME)
            else:
                xmls = sorted(p.glob("*.xml"))
                if not xmls:
                    raise FileNotFoundError(f"{p}: neither {SETTINGS_NAME} nor *.xml files")
                for x in xmls:
                    push(x, settings_override or _find_settings(p), None)
        elif p.suffix == ".settings":
            from_settings(p)
        elif p.suffix == ".source":
            from_source(p)
        elif p.suffix == ".population":
            for src in _read_population(p):
                from_source(src)
        elif p.suffix == ".xml":
            push(p, settings_override or _find_settings(p.parent), None)
        else:
            raise FileNotFoundError(f"{tok}: not a project root, .settings, .xml, .source or .population")
    return out


def _find_settings(processing_dir: Path) -> Optional[Path]:
    """``<root>/warp_tiltseries.settings`` for ``<root>/warp_tiltseries/x.xml``; else the single ``*.settings``
    in the parent."""
    cand = processing_dir.parent / SETTINGS_NAME
    if cand.exists():
        return cand
    hits = sorted(processing_dir.parent.glob("*.settings"))
    return hits[0] if len(hits) == 1 else None


def discover_series(
    xml_path: Path,
    settings_path: Optional[Path],
    source_params: Optional[dict],
    *,
    pix: Optional[float] = None,
    image_px: Optional[Tuple[int, int]] = None,
    volume_px: Optional[Tuple[int, int, int]] = None,
) -> WarpSeries:
    stem = xml_path.stem
    s = WarpSeries(stem=stem, xml_path=xml_path, settings_path=settings_path)

    settings = WarpSettings.from_file(str(settings_path)) if settings_path and settings_path.exists() else None
    s.settings = settings
    acq_pix = None
    if settings is not None:
        acq_pix = settings.pixel_size_a
        s.add("pixel_size_acquisition_a", acq_pix, f"{settings_path.name}#Import/PixelSize")
        bin_times = float(settings.get("Import", "BinTimes") or 0.0)
        if acq_pix:
            s.add("pix", acq_pix * (2.0**bin_times), f"{settings_path.name}#PixelSize x 2^BinTimes")
        if settings.tomo_dims_px and acq_pix:
            s.add("volume_dims_a", [d * acq_pix for d in settings.tomo_dims_px], f"{settings_path.name}#Tomo/Dimensions x PixelSize")
            s.add("volume_dims_px_acq", tuple(settings.tomo_dims_px), f"{settings_path.name}#Tomo/Dimensions")
        s.add("voltage", settings.voltage_kv, f"{settings_path.name}#CTF/Voltage")
        s.add("cs", settings.cs_mm, f"{settings_path.name}#CTF/Cs")
        s.add("amp_contrast", settings.amplitude_contrast, f"{settings_path.name}#CTF/Amplitude")
        s.add("exposure_per_tilt", settings.exposure_per_tilt, f"{settings_path.name}#Import/DosePerAngstromFrame")
        data_folder = settings.data_folder or "tomostar"
        s.tomostar_dir = (settings_path.parent / data_folder).resolve()
    if source_params:
        p = source_params
        if "PixelSize" in p and acq_pix is None:
            acq_pix = float(p["PixelSize"])
            s.add("pixel_size_acquisition_a", acq_pix, ".source#PixelSize")
            s.add("pix", acq_pix, ".source#PixelSize (BinTimes unknown)")
        if all(k in p for k in ("DimensionsX", "DimensionsY", "DimensionsZ")) and acq_pix:
            s.add("volume_dims_a", [float(p[k]) * acq_pix for k in ("DimensionsX", "DimensionsY", "DimensionsZ")], ".source#Dimensions x PixelSize")
        if "DosePerAngstromFrame" in p:
            s.add("exposure_per_tilt", -float(p["DosePerAngstromFrame"]), ".source#DosePerAngstromFrame")
    if s.tomostar_dir is None:
        # <root>/tomostar next to the processing folder, else the XML's DataDirectory attribute (often a
        # path from another machine), else the XML's own folder
        root = ElementTree.parse(xml_path).getroot()
        dd = root.get("DataDirectory")
        if (xml_path.parent.parent / "tomostar").is_dir():
            s.tomostar_dir = (xml_path.parent.parent / "tomostar").resolve()
        elif dd and Path(dd).is_dir():
            s.tomostar_dir = Path(dd).resolve()
        else:
            s.tomostar_dir = xml_path.parent.resolve()

    # explicit overrides for the XML reader
    if volume_px is not None:
        if pix is None and acq_pix is None:
            raise ValueError("--volume-px needs the acquisition pixel size (settings or --pix)")
        s.add("volume_dims_a", [v * (acq_pix or pix) for v in volume_px], "--volume-px x pixel", override=True)

    tomostar_path = s.tomostar_dir / f"{stem}.tomostar" if s.tomostar_dir else None
    if tomostar_path and tomostar_path.exists():
        s.tomostar_path = tomostar_path
        s.tomostar = WarpTomostar.from_file(str(tomostar_path))

    # image dims + tilt-image pixel from the header of a referenced average (what Warp reads)
    warp_probe = WarpAlignment.from_file(xml_path, pixel_size_a=pix or s.get("pix") or 1.0, strict_dims=False)
    avg = None
    for e in warp_probe.entries:
        if e.movie_path:
            movie = (s.tomostar_dir / e.movie_path.replace("\\", "/")).resolve() if s.tomostar_dir else None
            if movie is not None:
                cand = movie.parent / "average" / (movie.stem + ".mrc")
                if cand.exists():
                    avg = cand
                    break
    if avg is not None:
        h = mrc_header(avg)
        s.add("image_dims_px", (h["nx"], h["ny"]), f"{avg.parent.name}/{avg.name}#header")
        if h["voxel"][0] > 0:
            s.add("pix", round(h["voxel"][0], 6), f"{avg.parent.name}/{avg.name}#header", override=True)
    if image_px is not None:
        s.add("image_dims_px", tuple(image_px), "--image-px", override=True)
    if pix is not None:
        s.add("pix", pix, "--pix", override=True)
    xml_img = warp_probe.image_dimensions_physical
    if all(v > 0 for v in xml_img):
        s.add("image_dims_a", list(xml_img), f"{xml_path.name}#ImageDimensionsAngstrom")
    elif s.get("image_dims_px") and s.get("pix"):
        s.add("image_dims_a", [n * s.get("pix") for n in s.get("image_dims_px")], f"{s.source('image_dims_px')} x pix")
    xml_vol = warp_probe.volume_dimensions_physical
    if all(v > 0 for v in xml_vol):
        s.add("volume_dims_a", list(xml_vol), f"{xml_path.name}#VolumeDimensionsAngstrom", override=True)
    s.add("ctf_pixel_size_a", warp_probe.ctf_pixel_size_a, f"{xml_path.name}#CTF/PixelSize")
    for key, name in (("voltage", "Voltage"), ("cs", "Cs"), ("amp_contrast", "Amplitude")):
        if warp_probe.ctf_params.get(name):
            s.add(key, float(warp_probe.ctf_params[name]), f"{xml_path.name}#CTF/{name}")

    # reconstructions
    recon_dir = xml_path.parent / "reconstruction"
    if recon_dir.is_dir():
        s.reconstructions = sorted(p for p in recon_dir.glob(f"{stem}_*Apx.mrc"))
    return s
