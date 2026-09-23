# Astronomical Finding Chart Generator

A command-line Python tool for downloading public survey images and producing annotated astronomical finding charts with WCS coordinates, source markers, optional slit overlays, controlled sky orientation, and configurable image scaling.

The script supports **Pan-STARRS1**, **DESI Legacy Surveys**, **SkyMapper DR4**, **SDSS DR9**, **GALEX GR6/7**, **DES DR2**, **VISTA VIKING**, **VISTA Hemisphere Survey (VHS)**, **DSS2**, **2MASS**, and **AllWISE**.

Current release: **v1.9.7**.

### New in v1.9.7

- Reordered the `--survey` help choices so the compact survey names are shown first, followed by the long/canonical names.
- Added the short aliases `sm` = `skymapper`, `vik` = `viking`, and `aw` = `allwise`.


### New in v1.9.6

- Added **VISTA Hemisphere Survey (VHS)** as survey option `vhs`.
- VHS images are retrieved as FITS cutouts from the **ESO Science Archive** using TAP discovery plus the SODA cutout service.
- Supported VHS bands are `Y J H K`; `Ks` and `K_s` are accepted as aliases for `K`. The default VHS band is `K`.
- VHS SODA cutouts are native archive cutouts rather than HiPS2FITS visualization products, so they may be retained with `--fits`.

### New in v1.9.5

- The plotted sky field now remains fixed to the requested square FOV even when a long slit extends beyond the image boundary.
- Slit overlays are clipped at the finding-chart boundary and no longer contribute to Matplotlib autoscaling.

### New in v1.9.4

- The relative-geometry report now gives both the position angle from T to each secondary source and the equivalent opposite direction, 180 deg away.
- Whenever at least one secondary source is supplied, the finding-chart figure includes a compact box with the projected `DeltaRA` and `DeltaDec` offsets for **s2 - T**. This box is shown even when additional sources (`s3`, `s4`, ...) are present. The same arcsec/arcmin formatting rule is used as in the terminal report.

### New in v1.9.2

- Reformatted the relative-geometry terminal report for readability, with a blank line before the section and a separate multi-line block for each secondary source.

### New in v1.9.1

- When secondary sources are supplied, the program now reports their angular separation from the main target **T**, projected RA and Dec offsets, and the position angle from T to each source.
- Separations and coordinate offsets are shown in arcseconds, automatically switching to arcminutes when the absolute value exceeds 2 arcmin.
- Position angles are measured **North through East**, matching the convention used by `--angle`.

### New in v1.9.0

- Added **SDSS DR9** (`sdss`), **GALEX GR6/7** (`galex`), **DES DR2** (`des`), and **VISTA VIKING** (`viking`) through CDS HiPS2FITS.
- Replaced the previous **WISE All-Sky** option with **AllWISE** (`allwise`); the old `wise` survey option has been removed.
- No additional Python dependencies are required.

## Features

- Retrieve survey imaging from public services. Pan-STARRS1 and 2MASS prefer their standard backends and automatically fall back to CDS HiPS2FITS when needed; SDSS, GALEX, DES, VIKING, and AllWISE are retrieved directly through CDS HiPS2FITS; VHS is retrieved from the ESO Science Archive through TAP + SODA.
- Default field size of **4 x 4 arcmin**, with arbitrary square fields via `-f/--fov`.
- Accept sexagesimal or decimal RA/Dec coordinates.
- Plot up to **10 sources** in one chart.
- Per-source markers: `cross`, `circle`, `diamond`, `x`, or `slit`.
- Configurable marker colors and sizes.
- Configurable longslit length and width, with a circle marking the exact slit center.
- Use a controlled final orientation for every chart: PA = 0 deg (North up, East left) by default, or a requested position angle measured North through East.
- WCS RA/Dec axes, coordinate grid, and North/East compass.
- `auto`, `zscale`, and percentile intensity scaling.
- `asinh`, linear, square-root, and logarithmic stretches.
- Output to PNG, PDF, JPG, or EPS.
- Optionally keep downloaded FITS files from the normal survey backends with `--fits`; HiPS2FITS visualization products are intentionally not saved. Native ESO VHS SODA cutouts can be saved.
- Compact source information in the chart title, automatically limited to two lines.

## Requirements

Python **3.10 or newer** is required.

Install the required packages with:

```bash
python3 -m pip install numpy scipy matplotlib astropy astroquery requests
```

The program also requires internet access to retrieve survey images.

## Installation

Clone or download the repository, then make the script executable if desired:

```bash
chmod +x get_finding_chart.py
```

You can then run it as either:

```bash
python3 get_finding_chart.py [options]
```

or:

```bash
./get_finding_chart.py [options]
```

## Quick start

Only coordinates are required. The default is a **4 x 4 arcmin Pan-STARRS1 r-band** chart, displayed with **North up and East left**:

```bash
python3 get_finding_chart.py -c 17:05:35.520 -23:27:21.60
```

A different field size can be requested with `-f`:

```bash
python3 get_finding_chart.py \
    -c 17:05:35.520 -23:27:21.60 \
    -f 6
```

Run:

```bash
python3 get_finding_chart.py --help
```

for the complete command-line help.

## Coordinates and multiple sources

Coordinates may be given in sexagesimal form:

```bash
-c 17:05:35.520 -23:27:21.60
```

or decimal degrees:

```bash
-c 256.398 -23.456
```

Additional RA/Dec pairs can be appended to the same `-c` option:

```bash
python3 get_finding_chart.py \
    -c 17:05:35.520 -23:27:21.60 \
       17:05:36.100 -23:27:10.00 \
       17:05:34.800 -23:27:40.00
```

The first position is labelled **T**. Additional positions are labelled **s2**, **s3**, and so on. Up to 10 sources are supported.

Negative declinations can be supplied normally; they do not need special quoting.

## Surveys and bands

| Survey option | Survey | Supported bands | Default |
|---|---|---|---|
| `ps1` | Pan-STARRS1 | `g r i z y` | `r` |
| `legacy` | DESI Legacy Survey DR10 | `g r i z` | `r` |
| `skymapper` | SkyMapper DR4 | `u v g r i z` | `r` |
| `sdss` | SDSS DR9 | `u g r i z` | `r` |
| `galex` | GALEX GR6/7 | `FUV NUV` | `NUV` |
| `des` | DES DR2 | `g r i z Y` | `r` |
| `viking` | VISTA VIKING | `z Y J H K` | `K` |
| `vhs` | VISTA Hemisphere Survey (VHS) | `Y J H K` (`Ks`, `K_s` aliases) | `K` |
| `dss2` | DSS2 | `red/r`, `blue/b`, `ir/i` | `red` |
| `2mass` | 2MASS | `J H K` | `J` |
| `allwise` | AllWISE | `w1 w2 w3 w4` | `w1` |

Short survey aliases are `ls` = `legacy`, `sm` = `skymapper`, `vik` = `viking`, and `aw` = `allwise`. `panstarrs` is also accepted as an alias for `ps1`.

AllWISE also accepts the aliases `1`, `2`, `3`, `4`, or the approximate wavelengths `3.4`, `4.6`, `12`, and `22`. VIKING accepts `Ks`/`K_s` as aliases for `K`.

Example:

```bash
python3 get_finding_chart.py \
    -c 256.398 -23.456 \
    -s allwise -b w1 -f 8

python3 get_finding_chart.py \
    -c 256.398 -23.456 \
    -s vhs -b K -f 6
```

For the Legacy Survey, a different viewer layer can be selected with `--legacy-layer`.

For VHS, `K`, `Ks`, and `K_s` all select the ESO `Ks` filter. Availability of `Y` and `H` depends on the VHS sub-survey footprint; `J` and `Ks` provide the broadest VHS coverage.

> **SkyMapper note:** the SkyMapper service restricts individual cutouts to less than 10 arcmin on a side. Because every chart is reprojected to a controlled orientation, the input square is enlarged by a factor of `sqrt(2)` to avoid blank corners. The requested SkyMapper finding-chart field must therefore be smaller than about **7.07 arcmin** on a side.

## Image retrieval backends

The program uses different retrieval services for different surveys. Pan-STARRS1 and 2MASS keep the **preferred-backend first, HiPS2FITS only when needed** strategy. SDSS, GALEX, DES, VIKING, and AllWISE use CDS HiPS2FITS directly, while VHS uses the ESO Science Archive TAP + SODA services.

- **Pan-STARRS1:** the STScI Pan-STARRS cutout service is tried first. The returned FITS is checked over the field that will actually appear in the finding chart. If the preferred image fails to download or contains a substantial edge-connected/central region of non-finite pixels, the program retries with CDS **HiPS2FITS** using the Pan-STARRS DR1 HiPS for the requested band.
- **2MASS:** NASA SkyView is tried first. The final displayed field is checked in the same way; CDS **HiPS2FITS** is used only when the SkyView image fails or has inadequate coverage.
- **SDSS DR9, GALEX GR6/7, DES DR2, VISTA VIKING, and AllWISE:** CDS **HiPS2FITS** is the primary retrieval backend.
- **Legacy Survey:** the DESI Legacy Survey FITS cutout service is used directly.
- **DSS2:** NASA SkyView is used.
- **SkyMapper:** the SkyMapper DR4 SIAP service is used.
- **VHS:** the ESO Science Archive is queried through TAP and the selected public product is clipped with the ESO SODA service.

### Automatic coverage check for Pan-STARRS1 and 2MASS

The coverage test is deliberately applied to the **final chart field**, not simply to the larger downloaded input image. This avoids switching to HiPS merely because an enlarged cutout contains NaNs outside the region that will be displayed.

The current fallback criteria are conservative:

- more than **1%** of the final field consists of non-finite pixels connected to an image border; or
- more than **50%** of the small central test region is non-finite.

Isolated internal NaNs do not trigger the fallback. The program prints a diagnostic such as:

```text
STScI Pan-STARRS final-field coverage: 99.73% finite; edge-connected missing=0.27%; center missing=0.00%
```

If the preferred image is incomplete and HiPS2FITS also fails, the program keeps the usable-but-incomplete preferred image and prints a warning rather than failing the entire finding-chart request.

HiPS2FITS is useful near boundaries between individual survey images because the HiPS representation combines survey imaging into a hierarchical sky map before generating the requested cutout. These outputs are resampled visualization/cutout products and should not be treated as the preferred source for calibrated survey photometry when standard survey products are available.

## Relative source geometry

When more than one source is supplied with `-c/--coordinates`, the program prints a compact geometry summary near the end of the run, after any retrieval/output warnings. The first coordinate is the main target **T**; all additional sources are measured relative to it.

For each secondary source the report gives:

- total angular separation from T;
- projected on-sky RA offset (`DeltaRA`, positive toward East);
- Dec offset (`DeltaDec`, positive toward North);
- position angle `PA(T->sN)` in degrees, measured North through East;
- the equivalent opposite PA, 180 deg away, useful because a slit can be described in either direction along the same axis.

Angular values are reported in arcseconds unless their absolute value is greater than 2 arcmin, in which case arcminutes are used. `DeltaRA` is the projected on-sky longitudinal offset rather than the raw numerical difference in RA coordinates. The reported PA can be used directly with `--angle` to orient the chart/slit along the line joining T and the secondary source. The opposite PA is `(PA + 180) mod 360` and describes the same slit axis in the reverse direction. Whenever at least one secondary source is supplied, the `DeltaRA` and `DeltaDec` offsets for **s2 - T** are also shown in a small box inside the finding-chart figure, even if additional sources are present.

Example terminal output:

```text
Relative geometry from T:
  (Delta RA is the projected on-sky offset; +RA=east, +Dec=north.)

s2:
    separation=37.42 arcsec;
    DeltaRA=+31.08 arcsec; DeltaDec=-20.85 arcsec;
    PA(T->s2)=123.86 deg; opposite PA=303.86 deg

s3:
    separation=1.35 arcmin;
    DeltaRA=-42.70 arcsec; DeltaDec=+1.15 arcmin;
    PA(T->s3)=328.19 deg; opposite PA=148.19 deg
```

## Markers

Choose markers with `-m/--marker`:

```bash
-m cross slit circle
```

The marker list follows the same order as the coordinate pairs. For example:

```bash
python3 get_finding_chart.py \
    -c 17:05:35.520 -23:27:21.60 \
       17:05:36.100 -23:27:10.00 \
       17:05:34.800 -23:27:40.00 \
    -m cross slit circle
```

means:

- `T` -> cross
- `s2` -> slit
- `s3` -> circle

If fewer markers than sources are given, the remaining sources default to `cross`.

Available markers are:

```text
cross  circle  diamond  x  slit
```

### Marker colors

Use one color for all sources:

```bash
--color orange
```

or give one color per source:

```bash
--color red blue orange
```

Matplotlib color names and color specifications are accepted. Shell-sensitive values such as hexadecimal colors should be quoted, for example:

```bash
--color '#ff0000' '#00a0ff'
```

Without `--color`, the main target is red and secondary sources are blue.

### Marker sizes

Use one multiplicative size factor for all sources:

```bash
--marker-size 1.3
```

or specify factors source by source:

```bash
--marker-size 1.0 0.8 1.4
```

A factor of `1.0` is the default size. For a `slit` marker, `--marker-size` changes only the circle marking the slit center; it does not change the physical slit dimensions.

## Slit marker

A slit can be assigned to any source with:

```bash
-m slit
```

or, for example, only to the second source:

```bash
-m cross slit
```

The default slit dimensions are:

- length: **240 arcsec** (4 arcmin)
- width: **1 arcsec**

Set them independently with:

```bash
--slit-length 180 --slit-width 0.8
```

`--slit-length` is given in arcseconds and accepts any positive finite value. `--slit-width` is also in arcseconds and accepts values from **0.1 to 60 arcsec**.

The slit is vertical in the final chart. Because every chart is reprojected to a controlled sky orientation, a slit without `--angle` is aligned North-South (PA = 0 deg). When `--angle` is specified, the requested PA points upward and the slit is aligned with that PA.

All slit markers in a given chart use the same `--slit-length` and `--slit-width` values.

## Position angle and default orientation

Every finding chart is reprojected onto a controlled final WCS so its orientation does not depend on the survey retrieval backend.

If `--angle` is omitted, the display uses **PA = 0 deg**, which means:

- North is up;
- East is left.

Use `-a/--angle` to choose a different orientation:

```bash
-a 45
```

PA is measured in degrees **North through East**. The output image is reprojected so the requested PA points upward. The PA is shown in the title and encoded in the output filename only when `--angle` is explicitly supplied; the default North-up orientation does not add `PA0` to the filename.

Example:

```bash
python3 get_finding_chart.py \
    -c 17:05:35.520 -23:27:21.60 \
    -a 45 \
    -m slit \
    --slit-length 180 \
    --slit-width 0.8
```

## Image scaling and contrast

The default display uses:

```text
--scale auto
--stretch asinh
--low 1
--high 99.5
```

Available intensity-limit methods are:

- `auto` - sigma-clipped statistics constrained by `--low` and `--high` percentiles
- `zscale` - Astropy ZScale
- `percentile` - direct percentile limits from `--low` and `--high`

Available stretches are:

```text
asinh  linear  sqrt  log
```

For example:

```bash
python3 get_finding_chart.py \
    -c 17:05:35.520 -23:27:21.60 \
    --scale percentile --low 0.5 --high 99.8 \
    --stretch sqrt
```

Absolute intensity limits can be forced with `--vmin` and/or `--vmax`.

## Saving the downloaded FITS image

By default, the internal survey FITS file is temporary and is removed after the chart is produced.

Use:

```bash
--fits field.fits
```

to keep the downloaded FITS when the final chart is based on one of the normal survey/retrieval backends. For example:

```bash
python3 get_finding_chart.py \
    -c 17:05:35.520 -23:27:21.60 \
    -s legacy \
    --fits field.fits
```

For Pan-STARRS1 and 2MASS, the program first downloads and evaluates the preferred STScI/SkyView product in a temporary location. Only after that decision is made can a requested persistent FITS be written. This prevents an inadequate preferred image from being left behind when the chart subsequently switches to HiPS2FITS.

If the final chart uses **CDS HiPS2FITS**, the requested `--fits` file is deliberately **not created or overwritten**. This includes the SDSS, GALEX, DES, VIKING, and AllWISE options, which use HiPS2FITS as their primary backend. The program prints a warning explaining that the HiPS2FITS image is a finding-chart / visualization product rather than the preferred image for survey photometry, and suppresses persistent FITS output to reduce the risk of misuse.

If the preferred Pan-STARRS STScI or 2MASS SkyView image is retained, `--fits` is honored. If that preferred image is incomplete, HiPS2FITS is attempted, but HiPS2FITS itself fails, the program falls back to the incomplete preferred image; in that case `--fits` can still save that image and the warning in the terminal should be noted.

For **VHS**, `--fits` is also honored: the saved file is the native ESO SODA cutout returned by the archive, not a HiPS2FITS rendering.

The finding chart is always reprojected to its controlled final orientation. A saved FITS file is therefore the **original downloaded input cutout**, not the reprojected finding-chart image. The input is enlarged enough to support the final square reprojection without blank corners.

## Output formats and filenames

The default output format is PNG. Available formats are:

```text
png  pdf  jpg  eps
```

Select one with `-e/--extension`:

```bash
-e pdf
```

Use `-i/--id` to change the filename prefix:

```bash
-i GRB230307A
```

Output filenames include the identifier, survey, band, main-target coordinates, field size, and PA when applicable. A typical filename is:

```text
GRB230307A_ps1_r_256p398_-23p456_4arcmin_PA45.png
```

## Chart title and source information

The first title line contains general field information, for example:

```text
Pan-STARRS1 | band r | 4x4 arcmin | PA 45 deg
```

Source coordinates are packed into at most two additional information lines. Slit dimensions are written immediately after the coordinates of the source carrying the slit, for example:

```text
s2: 17:05:36.100 -23:27:10.00 [slit 180" x 0.8"]
```

The program fits as many complete source entries as possible. If there is not enough room to print every source coordinate, the text ends with:

```text
[...]
```

This truncation affects only the title: all supplied sources are still plotted on the finding chart.

## Complete example

```bash
python3 get_finding_chart.py \
    -c 17:05:35.520 -23:27:21.60 \
       17:05:36.100 -23:27:10.00 \
       17:05:34.800 -23:27:40.00 \
    -f 5 \
    -i target \
    -s ps1 -b r \
    -a 45 \
    -m cross slit circle \
    --color red blue orange \
    --marker-size 1.0 1.0 0.8 \
    --slit-length 180 \
    --slit-width 0.8
```

## Notes and limitations

- Survey availability and sky coverage are determined by the external survey services. HiPS mosaicking can remove individual-image boundary gaps, but it cannot create data outside the actual survey footprint. VHS availability also depends on the selected band and VHS sub-survey footprint.
- Pan-STARRS1 and 2MASS prefer STScI/SkyView respectively; HiPS2FITS is used automatically only when the preferred retrieval fails or the final displayed field has substantial missing edge/central coverage. SDSS, GALEX, DES, VIKING, and AllWISE use HiPS2FITS directly. VHS uses ESO TAP discovery and a native SODA FITS cutout.
- The automatic coverage thresholds are pragmatic finding-chart safeguards, not a scientific data-quality metric. A field passing the check is not thereby certified for photometry.
- If both the preferred retrieval and HiPS2FITS fail, the request fails. If the preferred image exists but is incomplete and the HiPS fallback fails, the program uses the incomplete preferred image and prints a warning.
- Every chart is reprojected, including the default North-up view, so interpolation is part of the finding-chart rendering process.
- The program does not redistribute survey data; images are retrieved on demand from the corresponding public services.
- HiPS2FITS products are intended here for finding-chart construction and visualization, not as replacements for native calibrated survey products used for photometry.
- When survey images are used in publications, presentations, or data products, cite or acknowledge the original survey according to its data-use policy.
- CDS requests that users of HiPS2FITS acknowledge the service as a tool provided by CDS, Strasbourg, France.
- Any persistent `--fits` image remains the original downloaded input cutout rather than the reprojected finding-chart image.

## Author

Andrea Rossi

Development assistance: OpenAI ChatGPT.
