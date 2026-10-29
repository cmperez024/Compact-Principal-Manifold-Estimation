# Reproduce the tuning-selection figures

## One-dimensional CV figure

Run `code/plot_1d_tuning_figures.R` directly. It generates its own missing
inputs and writes both periodic and interval CV figures beside this guide.

## Flower and cashew surface panels

Run `code/plot_sphere_panels.R` directly. Its embedded generator computes
both the flower and cashew persistence diagrams and PME fits. It creates its
own working directories and writes `flower_panel.png` and `cashew_panel.png`.
Neither plotter executes another simulation script. The separate generators
remain optional precomputation tools for these same experiments.

The production surface configuration is N=2500, noise SD 0.02 (variance
0.0004), and a 100-value lambda grid. The smaller N=250/10-lambda setup was
an older demonstration configuration. The archived production source and a
direct pixel comparison with the published point-cloud panel support N=2500
and SD 0.02; the latest draft's variance 0.004 is inconsistent with these.
Selection uses lambda strictly below the example's elbow cutoff, while all
available diagnostic values remain visible. Residual-average axes are labeled
"Volume-averaged conditional squared residual" as in the draft's TeX overlay.

The full calculations are expensive. Script paths can be supplied from any
working directory. For reduced validation set `PME_SMOKE_TEST=1` and direct
`PME_OUTPUT_DIR` to a temporary directory. `code/hpc/sphere.sh` preserves the
historical cluster settings and requires adapting its path before submission.
