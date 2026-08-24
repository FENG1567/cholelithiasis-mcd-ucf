# Data availability and provenance

## Ready-to-adapt manuscript statement

National Center for Health Statistics (NCHS) Multiple Cause of Death public-use mortality files for 1999–2024 and CDC WONDER query systems are publicly available from their respective official websites. The authors do not redistribute individual death records in this repository. The public GitHub repository <https://github.com/FENG1567/cholelithiasis-mcd-ucf> contains the analysis code, archived CDC WONDER query specifications, de-identified aggregate derived tables, figure source data, and release-integrity checks supporting the reported results. The current versioned release is v0.1.0 (<https://github.com/FENG1567/cholelithiasis-mcd-ucf/releases/tag/v0.1.0>) and is distributed under the MIT License. No DOI or public archive identifier has been assigned.

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

## Release status

- Public repository release: v0.1.0.
- Software license: MIT; see [`LICENSE`](../LICENSE).
- Repository citation metadata: [`CITATION.cff`](../CITATION.cff).
- No DOI or public archive identifier has been assigned.

No dataset DOI, embargo, or access committee is asserted here because none has
been supplied by the authors.
