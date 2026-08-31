# Reproduce the tuning-selection figures

## One-dimensional CV figure

Run `code/generate_1d_results.R`, then `code/plot_1d_tuning_figures.R`, with `results_1d/` and `job_figs_1d/` created under `code/`. The plotting script writes both periodic and interval CV figures.

## Flower and cashew surface panels

1. Create `code/results_sph/`, `code/panel_flower/`, and `code/panel_cashew/`.
2. Run `code/generate_sphere_results.R`.
3. Enable the historically commented flower TDA block so `results_sph/flower_tda.rds` is created. The cashew TDA block is active.
4. Run `code/plot_sphere_panels.R` to assemble the panels.

The sphere stage is HPC-scale. `code/hpc/sphere.sh` preserves the original SLURM script; adjust its `R CMD BATCH` path for the cluster location.
