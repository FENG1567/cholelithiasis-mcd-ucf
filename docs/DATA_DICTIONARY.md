# Data dictionary and file-to-result mapping

## Core variables

| Variable | Definition | Unit / values |
|---|---|---|
| `year` | Death year | 1999–2024 |
| `A_record_axis_K80` | Resident death records with record-axis K80 mention | count |
| `B_main_A_and_UCD_K80` | A records also having UCD K80 | count |
| `gap_A_minus_B` | A minus B | count |
| `UCF` | B divided by A | probability |
| `std_UCF` | Fixed-weight direct age-sex standardized UCF | probability |
| `CI_lo`, `CI_hi` | 95% confidence interval bounds | probability or contrast scale |
| `std_probability` | Standardized probability of a mutually exclusive UCD destination, conditional on gap | probability; annual total = 1 |

## Main figures and source data

| Result | Source data |
|---|---|
| Figure 1: study universe and estimands | `data/figure_source_data/Figure1_source_data.csv` |
| Figure 2: annual K80 UCF | `data/figure_source_data/Figure2_source_data.csv` |
| Figure 3: composition/decomposition | `data/figure_source_data/Figure3_source_data.csv` |
| Figure 4: recorded UCD destinations | `data/figure_source_data/Figure4_source_data.csv` |
| Figure 5: internal cross-disease benchmark | `data/figure_source_data/Figure5_source_data.csv` |
| Figure 6: robustness boundaries | `data/figure_source_data/Figure6_source_data.csv` |
| Supplementary Figures S1–S3 | corresponding `Supplementary_Figure_S*_source_data.csv` |

## Main tables

- `data/derived/Table1_Cohort_and_estimands.csv`: cohort boundary and frozen estimands.
- `data/derived/Table2_Main_estimates.csv`: central standardized estimates and endpoint contrasts.

## Important supplementary tables

- `k80_annual_main.csv`: annual A, B, gap, crude UCF and official-UCD reconciliation counts.
- `k80_standardized_annual.csv`, `k80_standardized_contrasts.csv`, `k80_standardization_weights.csv`: standardization outputs and weights.
- `ucd_destination_annual.csv`, `ucd_destination_ontology.csv`: mutually exclusive recorded-UCD allocation categories and ontology.
- `k80_composition_*.csv`, `k80_kitagawa_*.csv`: compositional and decomposition analyses.
- `cross_disease_*.csv`: exploratory internal comparison outputs.
- `robustness_matrix_current.csv`: executed robustness specifications.

All CSVs are UTF-8 encoded. Counts are unweighted record counts unless specified otherwise. No missing-value code is used in the released aggregate data; non-estimable quantities are represented by empty fields and described in table notes.
