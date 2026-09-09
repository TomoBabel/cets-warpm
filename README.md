# cets-warpm

Converters between Warp 2.0 / M tilt-series projects and the CETS cryo-ET standard, rigid profile
`cets-rigid/0.1` (see `cryoet-alignment/docs/cets.md`). M's refined geometry lives in the same
`warp_tiltseries/<stem>.xml` files; a `.source` or `.population` is accepted as a series enumerator.

```
cets-warpm to-cets warp/ -o out/warp.cets.json                # project root (warp_tiltseries.settings)
cets-warpm to-cets population/project.population -o out/m.cets.json
cets-warpm from-cets out/warp.cets.json -o warp_out/ --voltage 300 --cs 2.7 --amp-contrast 0.07 --dose-per-tilt 3.87
```

`to-cets` reads per series: the settings (acquisition pixel, `BinTimes`, `Tomo/Dimensions`, CTF constants,
`DataFolder`), the tomostar, the header of a referenced frame average (tilt-image size and pixel), the XML
(`Angles`, `AxisAngle`, `AxisOffsetX/Y`, `LevelAngleX/Y`, `UseTilt`, `Dose`, `MoviePath`, CTF grids) and the
reconstructions in `warp_tiltseries/reconstruction/`. Spatially constant deformation grids are folded into
the rigid shifts; spatially varying grids are refused unless `--drop-locals`. `AreAnglesInverted` and the
pixel sizes (acquisition / tilt image / CTF fitting) go to the companion.

`from-cets` writes `warp_tiltseries.settings` (once per project), `tomostar/<stem>.tomostar` (one row per
tilt image, darks included) and `warp_tiltseries/<stem>.xml` (rigid, `UseTilt=False` for images without
alignment, per-tilt CTF grids). Per-image X rotations must agree (Warp has one `LevelAngleX`) — the written
fields are reloaded and the rotation operators compared. Frames are not written; the hint names the
MoviePath targets.

| option | default / derivation |
|---|---|
| `--pix` | frame-average header › settings `PixelSize × 2^BinTimes`; required otherwise |
| `--image-px WxH`, `--volume-px XxYxZ` | XML attributes › average header / settings `Tomo/Dimensions`; required for zeroed XMLs |
| `--settings FILE` | `<root>/warp_tiltseries.settings` next to the processing folder |
| `--drop-locals`, `--no-ctf`, `--paths` | off / off / relative (warned) |
| from-cets `--voltage --cs --amp-contrast --dose-per-tilt --pix` | companion values; required otherwise (a settings file needs them) |
| from-cets `--frames-dir` | `ROOT/frames` (warned) |
| from-cets `--alignment NAME\|N`, `--tomogram ID` | required when a region has several |
