# Astronomical Finding Chart Generator

Command-line Python tool for downloading public survey images and creating annotated astronomical finding charts with WCS coordinates, source markers, optional slit overlays, controlled position angle, configurable image scaling, and Gaia-assisted acquisition for faint targets.

Current release: **v2.0.4**.

## Author

**Andrea Rossi**  
Development assistance: OpenAI ChatGPT.

## License

Released under the **MIT License**. See [LICENSE](LICENSE).

## Features

- Retrieve imaging from **Pan-STARRS1, DESI Legacy Survey DR10, SkyMapper DR4, SDSS DR9, GALEX GR6/7, DES DR2, VISTA VIKING, VISTA Hemisphere Survey (VHS), DSS2, 2MASS, and AllWISE**.
- Accept decimal or sexagesimal RA/Dec coordinates and plot up to **10 sources**.
- Default **4 x 4 arcmin** field, with arbitrary square FOVs through `-f/--fov`.
- Controlled chart orientation: North up / East left by default, or any PA with `-a/--angle`.
- Source markers: `x`, `cross`, `circle`, `diamond`, and `slit`.
- Default marker colors use a **color-blind-friendly Okabe-Ito palette**; `--color` overrides the defaults.
- North/East compass is always drawn in **dark green**, independently of the source-marker palette.
- Configurable slit length, width, marker size, and per-source colors.
- Long slits are clipped to the requested square FOV and do not change the plot size.
- Relative geometry for secondary sources: separation, projected `DeltaRA`, `DeltaDec`, PA, and opposite PA.
- Compact `s2 - T` offset box in normal mode.
- **Faint-source acquisition mode** using Gaia DR3:
  - `-F`, `--faint`, or `--faint-source`;
  - searches within **90 arcsec** of T;
  - default selection `G_RP < 19`;
  - Gaia positions propagated to epoch **2026.0** by default;
  - selected Gaia star becomes the chart center;
  - T defaults to a slit marker unless explicitly overridden;
  - automatic PA aligns the selected Gaia acquisition star with T;
  - `--faint N` selects the Nth closest suitable star;
  - optional custom `G_RP` range with `--faint MAG_MIN MAG_MAX` or `--faint N MAG_MIN MAG_MAX`;
  - `--epoch YEAR` changes the proper-motion epoch;
  - terminal output lists up to the five closest Gaia candidates, plus the selected star when needed;
  - only the selected acquisition star is plotted, with a `GN - T` offset box.
- Intensity scaling with `auto`, `zscale`, or percentile limits and `asinh`, linear, square-root, or logarithmic stretches.
- `--invert` for a dark-background grayscale display.
- Output to PNG, PDF, JPG, or EPS.
- Optional saving of suitable downloaded FITS cutouts with `--fits`.

## Requirements

Python **3.10 or newer**.

Install dependencies with:

```bash
python3 -m pip install -r requirements.txt
```

The program requires internet access to query survey services and Gaia DR3.

## Installation

Clone or download the repository. Optionally make the script executable:

```bash
chmod +x get_finding_chart.py
```

Run with:

```bash
python3 get_finding_chart.py [options]
```

or:

```bash
./get_finding_chart.py [options]
```

For all options:

```bash
python3 get_finding_chart.py --help
```

## Quick start

Default Pan-STARRS1 r-band chart:

```bash
python3 get_finding_chart.py -c 17:05:35.520 -23:27:21.60
```

Different field size and survey:

```bash
python3 get_finding_chart.py \
    -c 17:05:35.520 -23:27:21.60 \
    -f 6 -s legacy -b r
```

Multiple sources:

```bash
python3 get_finding_chart.py \
    -c 17:05:35.520 -23:27:21.60 \
       17:05:36.100 -23:27:10.00 \
       17:05:34.800 -23:27:40.00
```

The first source is labelled **T**; additional sources are `s2`, `s3`, etc.

## Surveys and bands

| Option | Survey | Bands | Default |
|---|---|---|---|
| `ps1` | Pan-STARRS1 | `g r i z y` | `r` |
| `legacy` / `ls` | DESI Legacy Survey DR10 | `g r i z` | `r` |
| `skymapper` / `sm` | SkyMapper DR4 | `u v g r i z` | `r` |
| `sdss` | SDSS DR9 | `u g r i z` | `r` |
| `galex` | GALEX GR6/7 | `FUV NUV` | `NUV` |
| `des` | DES DR2 | `g r i z Y` | `r` |
| `viking` / `vik` | VISTA VIKING | `z Y J H K` | `K` |
| `vhs` | VISTA Hemisphere Survey | `Y J H K` | `K` |
| `dss2` | DSS2 | `red/r blue/b ir/i` | `red` |
| `2mass` | 2MASS | `J H K` | `J` |
| `allwise` / `aw` | AllWISE | `w1 w2 w3 w4` | `w1` |

`panstarrs` is also accepted as an alias for `ps1`. VIKING and VHS accept `Ks` / `K_s` as aliases for `K` where applicable. AllWISE also accepts `1`, `2`, `3`, `4` and the approximate wavelengths `3.4`, `4.6`, `12`, `22`.

Pan-STARRS1 and 2MASS prefer their standard backends and can fall back to CDS HiPS2FITS when coverage is inadequate. SDSS, GALEX, DES, VIKING, and AllWISE use CDS HiPS2FITS directly. VHS is retrieved from the ESO Science Archive via TAP + SODA.

## Markers, colors, and slit

Choose one marker per source if desired:

```bash
-m x slit circle
```

Missing markers default to `x`.

Override default colors with one color for all sources:

```bash
--color orange
```

or one color per source:

```bash
--color '#D55E00' '#0072B2' '#009E73'
```

The default slit is **240 arcsec x 1 arcsec**. Change it with:

```bash
--slit-length 180 --slit-width 0.8
```

A user-specified `--marker` always overrides the automatic slit marker used for T in faint-source mode.

## Position angle and geometry

The default orientation is **PA = 0 deg**, i.e. North up and East left.

Set another PA with:

```bash
-a 45
```

PA is measured North through East. When multiple sources are supplied, the terminal reports their separation, projected RA/Dec offsets, PA, and opposite PA relative to T. Angular values are shown in arcsec and switch to arcmin above 120 arcsec.

## Faint-source acquisition mode

Use:

```bash
python3 get_finding_chart.py \
    -c 17:05:35.520 -23:27:21.60 \
    --faint
```

Default behavior:

- search Gaia DR3 within 90 arcsec of T;
- require `G_RP < 19`;
- propagate positions to epoch 2026.0;
- choose the closest suitable star as `G1`;
- center the chart on that star;
- align the automatic PA from the acquisition star toward T;
- plot T with a slit unless the user specified another marker.

Useful variants:

```bash
--faint 2
--faint 14 18
--faint 2 14 18
--epoch 2000
```

These select G2, restrict the `G_RP` range, combine both, or change the propagation epoch. An explicit `--angle` overrides the automatic Gaia-star PA.

The program warns, but does not modify the request, if the chosen FOV or slit is too short to include both T and the acquisition star.

## Image display

Defaults:

```text
--scale auto
--stretch asinh
--low 1
--high 99.5
```

Available scale modes:

```text
auto  zscale  percentile
```

Available stretches:

```text
asinh  linear  sqrt  log
```

Use `--invert` for a dark-background image. `--vmin` and `--vmax` can force absolute display limits.

## FITS and output files

Default output is PNG. Other formats:

```bash
-e pdf
-e jpg
-e eps
```

Use `-i/--id` to set the output filename prefix.

Use:

```bash
--fits field.fits
```

to retain a suitable downloaded FITS cutout. HiPS2FITS products are intentionally not saved as persistent FITS files because they are used here as finding-chart / visualization products rather than preferred calibrated survey images. Native ESO VHS SODA cutouts can be saved.

## Notes

- Survey availability and footprint are determined by the external services.
- Finding-chart images are reprojected to the requested orientation, so interpolation is part of the rendering process.
- Persistent FITS files, when allowed, are the downloaded input cutouts rather than the reprojected finding-chart image.
- Faint-source mode depends on Gaia DR3 astrometry and photometry; proper-motion propagation uses Gaia `pmra` and `pmdec` when available.
- When using survey data in publications, follow the acknowledgement and citation policy of the original survey/service.
