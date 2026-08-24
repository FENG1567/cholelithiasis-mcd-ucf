# Data availability and provenance

## Ready-to-adapt manuscript statement

National Center for Health Statistics (NCHS) Multiple Cause of Death public-use mortality files for 1999–2024 and CDC WONDER query systems are publicly available from their respective official websites. The authors do not redistribute individual death records in this repository. This repository contains the analysis code, archived CDC WONDER query specifications, de-identified aggregate derived tables, figure source data, and release-integrity checks supporting the reported results. A permanent repository release and DOI will be added by the authors before submission.

## Dataset-to-location map

| Supporting material | Access route | Location / action |
|---|---|---|
| NCHS Multiple Cause of Death public-use files (1999–2024) | Reused public source | Obtain from NCHS Vital Statistics Online; not redistributed here |
| CDC WONDER query definitions | Within repository | `queries/cdc_wonder_xml/` |
| Aggregate analytic tables | Within repository | `data/derived/` |
| Source data for main and supplementary figures | Within repository | `data/figure_source_data/` |
| Code and validation checks | Within repository | `code/` and `verify_release.py` |

## Provenance

The annual NCHS public-use records were streamed and reduced to aggregate age-by-sex, outcome, and destination counts. This package contains only those aggregate outputs, no person-level records. CDC WONDER XML files document the corresponding query configurations. Definitions, missingness handling, standardization weights, and ontology files are released in the derived tables and `docs/DATA_DICTIONARY.md`.

## Author actions before public release

- Create a versioned public repository release and archive it with a persistent DOI.
- Replace the placeholder repository URL and DOI in `CITATION.cff`.
- Choose a code licence only after confirming institutional/author rights.
- Ensure the repository record links the final article DOI and cites the NCHS/CDC sources.

No dataset DOI, repository identifier, licence, embargo, or access committee is asserted here because none has been supplied by the authors.
