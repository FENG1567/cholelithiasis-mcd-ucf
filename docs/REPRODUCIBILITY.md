# Reproducibility guide

## Release levels

This public package supports two distinct levels of checking.

1. **Released-result verification**: run `python verify_release.py`. This is fully executable from the repository root and checks the released aggregate tables, figure source data, XML queries, hashes, numerical identities, and release boundary.
2. **Full raw-data regeneration**: obtain NCHS annual Multiple Cause of Death public-use files directly from NCHS and use the archived CDC WONDER XML specifications. Full reruns require the official year-specific fixed-width layouts/crosswalk. Those external source files are not redistributed here. The analytic definitions are fixed in `code/analysis/ANALYSIS_SPECIFICATION.md`.

## Frozen primary definitions

- **A**: a resident death record with a record-axis ICD-10 condition beginning `K80`.
- **B**: `A=1` and underlying cause beginning `K80`.
- **UCF**: `B/A`; the arithmetic complement is `(A-B)/A`.
- **Primary standardization**: direct age-by-sex standardization using a fixed 2018–2024 A-distribution. Unknown age/sex do not enter standardization.
- **Years**: 1999–2024.

The release intentionally does not include a record-level mapping file, line numbers, or any table that could facilitate small-cell reconstruction.

## Recreate figures

Install R and the packages listed in `environment/R_environment.md`, then run:

```bash
Rscript code/figures/revision_figures.R .
```

The figure script accepts the repository root as its first argument. It reads only `data/derived/`; it writes 36 exports to `results/recreated_figures/`, nine recreated source-data CSVs to `results/recreated_figure_source_data/`, and export/numeric QA files to `results/qa/`. The main manuscript figures use the nine frozen source-data CSVs already released in `data/figure_source_data/`.

## Expected checks

`verify_release.py` verifies:

- 6 main figure and 3 supplementary figure source-data files;
- both manuscript primary tables;
- A=51,084 and B=27,514 across the annual K80 table;
- annual `gap=A-B` identity;
- nine mutually exclusive destination categories with annual standardized probabilities summing to 1;
- file hashes recorded in `docs/FILE_MANIFEST.csv`;
- no private-path, credential, raw-record, cache, or internal-review marker in released text files.

## Reuse limits

These outputs describe recorded cause-of-death data. They do not establish clinical diagnostic accuracy, preventability, treatment, or causal mechanisms. Consult the original manuscript for interpretation and limitations.
