# cets-warpm

Convert Warp 2 / M tilt-series projects to the CETS cryo-ET standard and back.

The package implements the rigid profile `cets-rigid/0.1` of CETS, documented in
[`cryoet-alignment/docs/cets.md`](https://github.com/uermel/cryoet-alignment/blob/uermel/cets/docs/cets.md).
It exchanges the rigid part of a Warp tilt-series model: per-tilt angles, tilt axis, image shifts, the
level angles and the per-tilt CTF. Deformation grids that vary in space are refused unless you ask to
drop them. Grids that are constant are rigid contributions and are folded into the shifts, never dropped.

- `cets-warpm to-cets` reads a Warp project, single XMLs, or an M `.source` / `.population`, and writes
  one CETS dataset.
- `cets-warpm from-cets` writes a Warp project (settings, tomostars, tilt-series XMLs) from a CETS dataset.

M refines the same `warp_tiltseries/<stem>.xml` files, so an M population is converted by enumerating
its sources; the geometry always comes from the XMLs.

## Installation

```bash
pip install git+https://github.com/TomoBabel/cets-warpm.git
```

The package depends on `cryoet-alignment >= 0.3.0` (the CETS codec) and on the pinned
`cets_data_model` commit listed in `pyproject.toml`. Until cryoet-alignment 0.3.0 is on PyPI, install it
from its branch first:

```bash
pip install git+https://github.com/uermel/cryoet-alignment.git@uermel/cets
```

## Expected source layout

### A Warp project

Point `to-cets` at the project root (the directory holding `warp_tiltseries.settings`). Every XML in
the processing folder named by the settings is one tilt series.

```
warp/                                    # project root
├── warp_tiltseries.settings             # required for a root: Import/PixelSize, BinTimes, Tomo/Dimensions,
│                                        #   CTF constants, DosePerAngstromFrame, DataFolder, ProcessingFolder
├── tomostar/                            # settings DataFolder
│   ├── TS_01.tomostar                   # movie names, tilt angles, axis angle, dose per tilt (row order)
│   └── TS_02.tomostar
├── warp_frameseries/
│   └── average/
│       ├── TS_01_001_0.00_....mrc       # frame averages: tilt-image size and pixel size from the header
│       └── ...                          #   of the first one a MoviePath points to
└── warp_tiltseries/                     # settings ProcessingFolder
    ├── TS_01.xml                        # required: Angles, AxisAngle, AxisOffsetX/Y, LevelAngleX/Y, UseTilt, Dose,
    │                                    #   MoviePath, CTF, GridCTF*, deformation grids (audited)
    ├── TS_02.xml
    └── reconstruction/
        ├── TS_01_10.00Apx.mrc           # optional: each <stem>_<voxel>Apx.mrc becomes a Tomogram entity
        └── TS_02_10.00Apx.mrc
```

`MoviePath` entries in the XML are relative to the tomostar folder (`../warp_frameseries/<movie>.eer`);
the average is looked up as `average/<movie>.mrc` next to the movie.

### An M population

```
population/
├── project.population                   # Sources/Source Path -> the .source files
└── ...
warp/
├── warp_tiltseries.settings
└── warp_tiltseries/
    ├── warp.source                      # the tomostars of this source; each maps to <stem>.xml in this folder
    ├── TS_01.xml                        # M writes its refined geometry back into these files
    └── TS_02.xml
```

Pass the `.population` (every source) or a single `.source`. The settings file is found next to the
processing folder, or given with `--settings`.

### Minimal input

A directory of XMLs, or single XML files, also work. The settings file is then searched one level up
(`<parent>/warp_tiltseries.settings`, or the single `*.settings` there). Whatever the XML, the settings
and the frame averages do not provide must be given: `--pix`, `--image-px`, `--volume-px`. XMLs re-saved
by Warp carry `ImageDimensionsAngstrom="0, 0"`, so the frame average header or `--image-px` is needed
for those.

## Warp to CETS

```bash
cets-warpm to-cets warp/ -o cets/warp.cets.json
```

What happens per series:

1. Settings, tomostar, XML and the first available frame average are read. Pixel sizes are kept
   distinct: acquisition (`Import/PixelSize`), tilt image (average header), CTF fitting (`CTF/PixelSize`).
2. Every XML row becomes a `TiltImage` (including `UseTilt=False` rows, which get no alignment).
   `nominal_tilt_angle` is the stage angle (Warp stores its negative).
3. Grids are audited. Constant `GridMovementX/Y` and constant `GridVolumeWarp` are folded into the
   per-tilt shifts. Spatially varying grids stop the conversion unless `--drop-locals`.
4. The rigid operator per tilt is `Rz(axis) · Ry(−angle − LevelAngleY) · Rx(LevelAngleX)`, expressed
   in the tilt image's and the reference volume's centred physical frames with Warp's half-pixel centre
   convention corrected on the way in.
5. The reference volume is Warp's own box (`Tomo/Dimensions` x acquisition pixel, or the XML's
   `VolumeDimensionsAngstrom`) as a virtual tomogram `<stem>_volume`. Each file in `reconstruction/`
   becomes a tomogram `<stem>_tomo_<voxel>` at its header voxel size.
6. `AreAnglesInverted`, the CTF pixel size, the dose and the movie names go to the companion.

Console output (abridged):

```
== Position_16_3
   pix = 1.54  [discovered]  (average/Position_16_3_031_-45.00_20240716_210415_EER.mrc#header)
   volume_dims_a = [6307.84, 6307.84, 3080.0]  [discovered]  (warp_tiltseries.settings#Tomo/Dimensions x PixelSize)
   paths = 'relative'  [default]
   dose_per_tilt = 3.87  [discovered]  (warp_tiltseries.settings#Import/DosePerAngstromFrame)
   voltage = 300.0  [discovered]  (warp_tiltseries.settings#CTF/Voltage)
   cs = 2.7  [discovered]  (warp_tiltseries.settings#CTF/Cs)
   amp_contrast = 0.07  [discovered]  (warp_tiltseries.settings#CTF/Amplitude)
   [ok ] rows value=31 expected=31
   WARNING: paths defaulted to 'relative'; set it with --paths or config key series.Position_16_3.paths
wrote cets/warp.cets.json (1 region(s)) + warp.cets-companion.json
report: cets/warp.cets.report.json
```

A series with local deformation is refused with the grid names, so you can decide:

```
== Position_16_3
   ERROR: Position_16_3: spatially varying deformation grids (GridMovementX, GridMovementY); the rigid profile cannot represent them - pass --drop-locals to keep the rigid part only
```

Output:

```
cets/
├── warp.cets.json             # one Region per series: tilt series, movie stacks, alignment "warp", tomograms
├── warp.cets-companion.json   # acquisition/CTF pixel sizes, kV/Cs/amp, dose, AreAnglesInverted, UseTilt,
│                              #   movie names, folded grid values, what was dropped
└── warp.cets.report.json      # per series: provenance of every value, gates, warnings, errors
```

More examples:

```bash
# an M population, every source
cets-warpm to-cets population/project.population -o cets/m.cets.json

# one M source with an explicit settings file
cets-warpm to-cets warp/warp_tiltseries/warp.source --settings warp/warp_tiltseries.settings -o cets/m.cets.json

# single XMLs from a project whose settings are elsewhere
cets-warpm to-cets warp/warp_tiltseries/TS_01.xml warp/warp_tiltseries/TS_02.xml \
    --settings warp/warp_tiltseries.settings -o cets/two.cets.json

# re-saved XMLs without dimensions and no frame averages at hand
cets-warpm to-cets xmls/ -o cets/x.cets.json --pix 1.54 --image-px 4096x4096 --volume-px 4096x4096x2000

# keep the rigid part of a deformed series
cets-warpm to-cets warp/ -o cets/warp.cets.json --drop-locals
```

## CETS to Warp

```bash
cets-warpm from-cets cets/warp.cets.json -o warp_out/
```

What happens per region:

1. The alignment and the reference tomogram are selected (flags required only when a region has
   several).
2. Representability is checked. Warp has one `LevelAngleX` per series, so every tilt's operator must
   factor as `Rz(axis_t) · Ry(−angle_t) · Rx(LevelAngleX)` with a common X angle. Otherwise the series is
   refused.
3. `warp_tiltseries.settings` is written once per project (acquisition pixel, `Tomo/Dimensions` from the
   reference box, dose, kV/Cs/amp). Every further series must agree with it, or it is refused with the
   list of conflicts.
4. The tomostar gets one row per tilt image, dark ones included. The XML gets `UseTilt=False` for images
   without an alignment, per-tilt CTF grids, and the rigid fields with the centre convention corrected on
   the way out.
5. The written XML is reloaded and its rotation operators compared with the document's.

Frames are not written. The `MoviePath` entries point at `--frames-dir` (default `ROOT/frames`), and the
hint names the first missing target.

Console output (abridged):

```
== Position_16_3
   reference_tomogram = 'Position_16_3_volume'  [discovered]  (companion)
   pix = 1.54  [companion]
   voltage = 300.0  [companion]
   cs = 2.7  [companion]
   amp_contrast = 0.07  [companion]
   dose_per_tilt = 3.87  [companion]  (Warp holds one DosePerAngstromFrame per project)
   angles_inverted = False  [companion]
   frames_dir = 'warp_out/frames'  [default]
   [ok ] operators_rebuilt value=3.2e-15 expected='<= 1e-06 deg'
   [ok ] settings  settings written
   [ok ] rows_match_tomostar value=(31, 31) expected=31
   WARNING: frames_dir defaulted to 'warp_out/frames'; set it with --frames-dir or config key series.Position_16_3.frames_dir
   > place the movies (or their average/<name>.mrc) at the MoviePath targets, e.g. warp_out/tomostar/../frames/Position_16_3_031_-45.00_20240716_210415_EER.eer
   > WarpTools ts_reconstruct --settings warp_out/warp_tiltseries.settings --angpix <voxel> --perdevice 1
report: warp_out/cets_warpm.report.json
```

Output:

```
warp_out/
├── warp_tiltseries.settings   # once per project
├── tomostar/
│   └── TS_01.tomostar         # one row per tilt image
├── warp_tiltseries/
│   └── TS_01.xml              # rigid model, UseTilt, CTF grids 1x1xT, other grids 1x1x1
└── cets_warpm.report.json     # provenance, gates, hints per series
```

More examples:

```bash
# a document written by cets-aretomo3: the companion has no Warp constants, supply them
cets-warpm from-cets cets/at3.cets.json -o warp_out/ --voltage 300 --cs 2.7 --amp-contrast 0.07 \
    --dose-per-tilt 3.0 --frames-dir /data/frames

# a portal document, one region, explicit alignment
cets-warpm from-cets cets/10445.cets.json -o warp_out/ --region TS_105_5 --alignment portal18924

# the same with a config file instead of flags
cets-warpm from-cets cets/at3.cets.json -o warp_out/ --config cets.yaml
```

## Command reference

### `cets-warpm to-cets`

```
cets-warpm to-cets [OPTIONS] SOURCES...
```

`SOURCES` are project roots, `.settings` files, `.xml` files or directories of them, M `.source` files or
`.population` files, in any mix. Duplicates are removed.

| Option | Meaning | Derived from, when absent |
|---|---|---|
| `-o, --output FILE` | CETS dataset JSON to write (required) | |
| `--name TEXT` | dataset name | the output stem |
| `--settings FILE` | `warp_tiltseries.settings` to use for every XML | `<root>/warp_tiltseries.settings` next to the processing folder |
| `--pix Å` | tilt-image (frame average) pixel size | average header, then settings `PixelSize x 2^BinTimes`; error otherwise |
| `--image-px WxH` | tilt image size | XML `ImageDimensionsAngstrom`, then the average header; error otherwise |
| `--volume-px XxYxZ` | reconstruction box in acquisition pixels | XML `VolumeDimensionsAngstrom`, then settings `Tomo/Dimensions`; error otherwise |
| `--drop-locals` | keep the rigid part of a deformed series | off: spatially varying grids are refused |
| `--no-ctf` | ignore the CTF values | off |
| `--paths relative\|absolute` | how file paths are written into the document | `relative`, with a warning |
| `--voltage kV`, `--cs mm`, `--amp-contrast F` | companion values only | settings / XML CTF constants |
| `--dose-per-tilt D` | per-image exposure (e/Å²), companion only | settings `DosePerAngstromFrame` |
| `--fail-fast` | stop at the first failing series | continue, exit 1 at the end |
| `--overwrite` | replace existing outputs | error when outputs exist |
| `--config FILE` | YAML config with overrides (see below) | |

### `cets-warpm from-cets`

```
cets-warpm from-cets [OPTIONS] DOCUMENT
```

`DOCUMENT` is a CETS dataset JSON. Its companion manifest is read when present.

| Option | Meaning | Derived from, when absent |
|---|---|---|
| `-o, --output DIR` | Warp project root (required) | |
| `--region ID` | region(s) to convert; repeatable | all |
| `--alignment NAME\|N` | alignment instance name or 0-based index | required when a region has several |
| `--tomogram ID` | reference tomogram id | companion; required when a region has several |
| `--pix Å` | acquisition pixel size (settings `Import/PixelSize`) | companion; the document's tilt-image pixel otherwise |
| `--frames-dir DIR` | where the movies will live (MoviePath targets) | `ROOT/frames`, with a warning |
| `--dose-per-tilt D` | exposure per tilt (settings `DosePerAngstromFrame`) | companion; error otherwise |
| `--voltage kV`, `--cs mm`, `--amp-contrast F` | settings CTF constants | companion; error otherwise |
| `--angles-inverted` | Warp `AreAnglesInverted` (defocus handedness) | companion; off otherwise |
| `--no-ctf` | do not write CTF grids | off |
| `--fail-fast`, `--overwrite`, `--config FILE` | as above | |

## Values, defaults and the config file

Every value a command needs resolves through one chain:

```
CLI flag  >  --config FILE  >  companion manifest  >  discovered (settings, headers, XML)  >  package default
```

Each resolution is printed as `option = value  [source]` and recorded in the report. Falling back to a
package default always prints a warning that names the flag and the config key. Values without a sane
default (kV, Cs, amplitude contrast and dose for a settings file; pixel size and boxes on import) are
errors that name the flag.

```yaml
cets:                          # every package and command
  voltage: 300
  cs: 2.7
  amp_contrast: 0.07
  paths: relative
cets-warpm:
  to-cets: {drop_locals: false}
  from-cets: {frames_dir: /data/frames, dose_per_tilt: 3.87, angles_inverted: false}
series:                        # per tilt series; wins over the command section
  TS_01: {pix: 1.54, image_px: 4096x4096, volume_px: 4096x4096x2000}
```

```bash
cets-warpm to-cets warp/ -o cets/warp.cets.json --config cets.yaml
```

## What is refused, what is dropped

- Spatially varying deformation grids, unless `--drop-locals`. Constant grids are folded, never dropped.
- XMLs with zero dimensions and no frame average, unless `--image-px` / `--volume-px` are given.
- On export, per-tilt X rotations that do not share one `LevelAngleX` (tolerance 1e-3 degrees).
- On export, a second series whose settings disagree with the project's settings file.
- `GridAngleX/Y/Z` (particle-orientation grids) are recorded in the companion, not converted.
- `AreAnglesInverted` is carried, never applied: it affects Warp's defocus channel, not the geometry.

## Development

```bash
pip install -e '.[dev]'
pytest
pre-commit run --all-files     # ruff 0.11.12, ruff-format, mypy 1.8.0
```

`tests/test_cross_tool.py` converts a real AreTomo3 `.aln` through CETS into a Warp XML and checks that
independent reference implementations of the Warp and AreTomo3 projection models project 3D points
identically. It skips when those reference models are not installed. The checks of constant-grid folding
and representability live in cryoet-alignment (`tests/test_cets_warp_goldens.py`).
