# Public analysis specification

This document describes the fixed analysis logic for raw-data regeneration. It deliberately avoids retaining local paths, individual identifiers, or private audit examples.

1. Use NCHS annual Multiple Cause of Death public-use files for 1999–2024 and retain resident records.
2. Identify K80 on the **record axis**. Define A as any such record; define B only when A is present and the UCD begins K80. Do not substitute all official UCD-K80 deaths for B.
3. Produce age-by-sex aggregate counts; merge 0–14 and 15–24 years to 0–24 before fixed-weight standardization. Exclude unknown age/sex from standardization but retain them in total counts and missingness audits.
4. Use the 2018–2024 A-distribution as the primary fixed age-by-sex weight. Report 1999 and whole-period weights as sensitivity analyses. Do not renormalize weights by year.
5. In gap records (`A-B`), classify the UCD into the released nine-category ontology (`data/derived/ucd_destination_ontology.csv`). Categories must be mutually exclusive, exhaustive, and sum to one within each year after standardization.
6. Interpret all findings as recorded attribution patterns. Do not infer clinical correctness, preventability, discrimination, treatment effects, or causal coding changes.

The released result tables and XML query specifications allow independent inspection. A full raw-data rerun additionally requires the official file layouts/crosswalk, which should be downloaded directly from NCHS and versioned by the reproducing researcher.
