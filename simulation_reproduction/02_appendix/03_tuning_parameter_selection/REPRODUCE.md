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

The full calculations are expensive. Script paths can be supplied from any
working directory. For reduced validation set `PME_SMOKE_TEST=1` and direct
`PME_OUTPUT_DIR` to a temporary directory. `code/hpc/sphere.sh` preserves the
historical cluster settings and requires adapting its path before submission.
