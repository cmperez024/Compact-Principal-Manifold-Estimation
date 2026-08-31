# Reproduce `interval_S1_full.png`

The local `code/` directory contains the complete two-stage 1D pipeline.

1. Run `code/generate_1d_results.R` to create `results_1d/*.rds`.
2. Run `code/plot_1d_figures.R` to create `job_figs_1d/interval_S1_full.png` and the related tuning plots.

From this directory in R:

```r
setwd("code")
dir.create("results_1d", showWarnings = FALSE)
dir.create("job_figs_1d", showWarnings = FALSE)
source("generate_1d_results.R")
source("plot_1d_figures.R")
```

Required packages are declared at the top of the two scripts. This can be a long computation.
