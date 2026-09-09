"""``cets-warpm`` command line: ``to-cets`` (Warp/M project -> CETS) and ``from-cets`` (CETS -> Warp project)."""

import sys
from pathlib import Path

import click
from cryoet_alignment.io.cets.cli_support import (
    Report,
    SeriesReport,
    common_options,
    echo,
    finish,
    load_config,
    make_resolver,
    parse_alignment_selector,
    parse_size,
    print_series,
    selection_options,
)
from cryoet_alignment.io.cets.companion import Companion
from cryoet_alignment.io.cets.entities import dataset_entity, dump_json, load_dataset, validate_document

from cets_warpm import __version__
from cets_warpm.discover import discover_series, enumerate_xmls
from cets_warpm.from_cets import cets_to_warp
from cets_warpm.to_cets import warp_to_cets

PACKAGE = "cets-warpm"
TO_CETS_OPTIONS = {
    "pix",
    "image_px",
    "volume_px",
    "drop_locals",
    "no_ctf",
    "paths",
    "voltage",
    "cs",
    "amp_contrast",
    "dose_per_tilt",
    "settings",
}
FROM_CETS_OPTIONS = {
    "pix",
    "frames_dir",
    "tomo_size",
    "dose_per_tilt",
    "voltage",
    "cs",
    "amp_contrast",
    "angles_inverted",
    "no_ctf",
}


@click.group()
@click.version_option(__version__)
def main():
    """Warp / M <-> CETS (rigid profile cets-rigid/0.1)."""


@main.command("to-cets")
@click.argument("sources", nargs=-1, required=True)
@click.option(
    "-o", "--output", "output", required=True, type=click.Path(dir_okay=False), help="CETS dataset JSON to write."
)
@click.option("--name", default=None, help="Dataset name (default: the output stem).")
@click.option(
    "--settings",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="warp_tiltseries.settings to use for every XML.",
)
@click.option("--pix", type=float, default=None, help="Tilt-image (frame average) pixel size Å/px.")
@click.option("--image-px", "image_px", default=None, help="Tilt image size WxH (re-saved XMLs carry zeros).")
@click.option("--volume-px", "volume_px", default=None, help="Reconstruction box XxYxZ in acquisition pixels.")
@click.option(
    "--drop-locals", "drop_locals", is_flag=True, default=None, help="Keep the rigid part of a deformed series."
)
@click.option("--no-ctf", "no_ctf", is_flag=True, default=None)
@click.option("--paths", type=click.Choice(["relative", "absolute"]), default=None)
@click.option("--voltage", type=float, default=None, help="kV (companion only).")
@click.option("--cs", type=float, default=None, help="mm (companion only).")
@click.option("--amp-contrast", "amp_contrast", type=float, default=None)
@click.option(
    "--dose-per-tilt", "dose_per_tilt", type=float, default=None, help="Per-image exposure e/Å² (companion only)."
)
@common_options
def to_cets(sources, output, name, config_path, overwrite, fail_fast, **cli):
    """Convert Warp tilt series (project root, .settings, .xml, M .source/.population) to a CETS dataset."""
    out = Path(output)
    if out.exists() and not overwrite:
        raise click.ClickException(f"{out} exists (use --overwrite)")
    out.parent.mkdir(parents=True, exist_ok=True)
    config = load_config(config_path, TO_CETS_OPTIONS)
    flags = {k: v for k, v in cli.items() if v is not None}
    image_px = parse_size(flags.pop("image_px", None), 2, "image-px")
    volume_px = parse_size(flags.pop("volume_px", None), 3, "volume-px")
    settings_override = Path(flags.pop("settings")) if "settings" in flags else None
    report = Report(PACKAGE, "to-cets")
    companion = Companion(generator=f"{PACKAGE} {__version__}")
    regions = []
    try:
        triples = enumerate_xmls(list(sources), settings_override)
    except FileNotFoundError as e:
        raise click.ClickException(str(e)) from e
    for xml, settings, src in triples:
        sr = SeriesReport(xml.stem)
        report.series.append(sr)
        if src and "__missing__" in src:
            sr.warnings.append(f"{src['__missing__']}: listed in the .source but {xml.name} does not exist - skipped")
            print_series(sr)
            continue
        try:
            series = discover_series(xml, settings, src, pix=flags.get("pix"), image_px=image_px, volume_px=volume_px)
            sr.warnings.extend(series.warnings)
            res = make_resolver(PACKAGE, "to-cets", flags, config, series.stem, sr)
            result = warp_to_cets(series, res, sr, out_dir=out.parent)
        except Exception as e:  # noqa: BLE001
            sr.error = str(e)
            print_series(sr)
            if fail_fast:
                break
            continue
        regions.append(result.region)
        companion.tilt_series[series.stem] = result.tilt_series_companion
        companion.alignments.append(result.alignment_companion)
        companion.tomograms.update(result.tomogram_companions)
        print_series(sr)
    if regions:
        ds = dataset_entity(name or out.name.split(".")[0], regions)
        validate_document(ds)
        dump_json(ds, out)
        companion.dump(Companion.path_for(out))
        echo(f"wrote {out} ({len(regions)} region(s)) + {Companion.path_for(out).name}")
    finish(report, out.with_name(out.name.split(".")[0] + ".cets.report.json"))


@main.command("from-cets")
@click.argument("document", type=click.Path(exists=True, dir_okay=False))
@click.option("-o", "--output", "output", required=True, type=click.Path(file_okay=False), help="Warp project root.")
@selection_options
@click.option("--pix", type=float, default=None, help="Acquisition pixel size Å/px (settings Import/PixelSize).")
@click.option(
    "--frames-dir",
    "frames_dir",
    type=click.Path(file_okay=False),
    default=None,
    help="Where the movies will live [ROOT/frames].",
)
@click.option(
    "--dose-per-tilt",
    "dose_per_tilt",
    type=float,
    default=None,
    help="Exposure per tilt e/Å² (settings DosePerAngstromFrame).",
)
@click.option("--voltage", type=float, default=None)
@click.option("--cs", type=float, default=None)
@click.option("--amp-contrast", "amp_contrast", type=float, default=None)
@click.option(
    "--angles-inverted",
    "angles_inverted",
    is_flag=True,
    default=None,
    help="Warp AreAnglesInverted (defocus handedness).",
)
@click.option("--no-ctf", "no_ctf", is_flag=True, default=None)
@common_options
def from_cets(document, output, regions, alignment, tomogram, config_path, overwrite, fail_fast, **cli):
    """Write a Warp project (settings, tomostar, rigid tilt-series XML) for the regions of a CETS dataset."""
    doc = Path(document)
    ds = load_dataset(doc)
    companion = Companion.load_for(doc)
    root = Path(output)
    config = load_config(config_path, FROM_CETS_OPTIONS)
    flags = {k: v for k, v in cli.items() if v is not None}
    report = Report(PACKAGE, "from-cets")
    wanted = set(regions)
    for region in ds.regions:
        if wanted and region.id not in wanted:
            continue
        sr = SeriesReport(region.id)
        report.series.append(sr)
        try:
            res = make_resolver(PACKAGE, "from-cets", flags, config, region.id, sr)
            cets_to_warp(
                region,
                res,
                sr,
                root=root,
                doc_dir=doc.parent,
                companion=companion,
                alignment_selector=parse_alignment_selector(alignment),
                tomogram_selector=tomogram,
                overwrite=overwrite,
            )
        except Exception as e:  # noqa: BLE001
            sr.error = str(e)
        print_series(sr)
        if sr.error and fail_fast:
            break
    root.mkdir(parents=True, exist_ok=True)
    finish(report, root / "cets_warpm.report.json")


if __name__ == "__main__":
    sys.exit(main())
