# R environment for figure recreation

Validated figure script dependencies:

- R >= 4.3
- ggplot2
- dplyr
- tidyr
- scales
- patchwork
- svglite
- ragg
- grid (base R)

Example installation:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "scales", "patchwork", "svglite", "ragg"))
```

The R script recreates graphics from released aggregate tables; it does not require or download NCHS individual-level records.
