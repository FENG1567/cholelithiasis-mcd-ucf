# Reproducibility package

This repository accompanies the study **Conditional Underlying-Cause Fractions Among U.S. Deaths With Record-Axis Cholelithiasis Codes, 1999–2024: Decomposition and Exploratory Internal Comparisons**.

It contains the aggregate analysis outputs needed to inspect the reported results, recreate publication figures, and run release-integrity checks. It does **not** contain death-certificate microdata, personal paths, credentials, internal review material, or record-level audit examples.

The repository is public and the current versioned release is [v0.1.0](https://github.com/FENG1567/cholelithiasis-mcd-ucf/releases/tag/v0.1.0), released under the MIT License. No DOI or public-archive identifier has been assigned. For a manuscript-ready code-availability statement, see [docs/CODE_AVAILABILITY.md](docs/CODE_AVAILABILITY.md).

## Quick start

Python 3.10 or newer is required. On Windows, use `py -3`; on macOS/Linux, use `python3`.

```bash
# Windows
py -3 verify_release.py

# macOS/Linux
python3 verify_release.py
```

The command verifies file hashes, expected source-data coverage, core numerical identities, destination closure, and prohibited-content patterns. It does not download or redistribute restricted/large source files.

To recreate the publication graphics after installing the documented R packages:

```bash
Rscript code/figures/revision_figures.R .
```

Generated graphics are written to `results/recreated_figures/`, recreated source-data CSVs to `results/recreated_figure_source_data/`, and checks to `results/qa/`. Their values should be compared with the frozen files in `data/figure_source_data/`.

## Package map

- `data/derived/`: de-identified aggregate result tables; no individual death records.
- `data/figure_source_data/`: source data for all 6 main and 3 supplementary figures.
- `queries/cdc_wonder_xml/`: archived CDC WONDER query specifications and manifest.
- `code/`: figure recreation and public validation code.
- `docs/`: definitions, provenance, data access, and reproducibility instructions.
- `environment/`: Python and R environment specifications.

## Data access and scope

NCHS Multiple Cause of Death public-use files and CDC WONDER are reused public sources. Obtain the original records directly from the official providers; do not rely on this repository as a source for individual-level mortality data. This package releases only derived aggregate values required to support the manuscript.

This release is distributed under the MIT License; see [`LICENSE`](LICENSE). The repository contains no individual-level mortality records, and no DOI or public-archive identifier has been assigned. See [docs/DATA_AVAILABILITY.md](docs/DATA_AVAILABILITY.md) for the manuscript-ready data-availability statement.
