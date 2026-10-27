# Changelog

## 5.1.0

- Add `update_github.sh` to copy a packaged release into a fresh clone of an existing Git repository, commit changed release files, and push its current branch.

## 5.0.2

- Keep the final target T visible in the chart's coordinate heading when many source blocks exceed its two-line layout.
- Include the README example image in the release archive.

## 5.0.1

- Center normal multi-source charts and downloaded fields on the first supplied source S, while keeping the final source T as the offset destination and output-name coordinate. Faint-source mode continues to center on the selected Gaia star.

## 5.0.0

- Interpret the final coordinate pair as T and preceding pairs as S, S2, S3, etc. Center normal charts and name outputs from the final target; calculate source offsets toward it.
- Apply the faint-source default slit and Gaia acquisition search to the final target. This is a major release because existing multi-coordinate commands change target unless their coordinate pairs are reordered.

## 4.0.1

- Correct input-source labels to `T`, `S`, `S2`, `S3`, ... consistently across chart markers, coordinate headings, offsets, terminal reports, and CLI help. Gaia acquisition labels are unchanged.

## 4.0.0

- Present accepted survey aliases as short/long pairs in `--help`.
- Add `-o/--output` for the chart; infer its format from the suffix and reject conflicts with `-e/--extension`.
- Add `--overwrite` for charts and requested native FITS files; existing outputs are now refused by default.
- Stage chart and retained FITS files beside their destinations before publication.
- This is a major release because repeated commands no longer replace existing charts or retained FITS files unless requested.

## 3.0.0

- Change the VHS default image band from K (ESO Ks) to J. Explicit `-b K` still requests Ks.
- This is a major release because an unchanged `-s vhs` command now produces a different band.

## 2.2.0

- Add UKIRT Hemisphere Survey DR1 J-band images via CDS HiPS2FITS (`-s uhs`).

## 2.1.2

- Add `-v` as a short form of `--version`.
- Package the missing 2.1.1 release and the current 2.1.2 release with documentation, dependencies, license, and tests.

## 2.1.1

- Brighten the auto-scaled sky background further by adjusting the lower cut while retaining the upper cut used for black star cores.
- Reduce chart raster resolution from 147 to 98 dpi.

## 2.1.0

- Lighten a dark sky background in default auto scaling using a quick sky-median estimate; preserve explicit absolute cuts and the other scale modes.
- Match the compass color to the Gaia acquisition marker blue.
- Reduce chart raster resolution from 220 to 147 dpi.
- Add `-g/--graphics` to open a successfully saved chart with the Linux or macOS default viewer.
