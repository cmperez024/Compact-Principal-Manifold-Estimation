# Reproduce the additional simulations

- `1d3D_sine_spiral.png` and `1d3D_star_moon.png`: run `code/one_dimensional_3d_examples.R`.
- `dist_lambda_infty.png`: run `code/asymptotic_simulations.R`.

Each script imports its own packages and embeds its own data generators.
Neither script reads another simulation's results. The 1D-in-3D script uses
the current CompactPME fitters and emits only the two advertised figures.
The asymptotic script computes the interval, circle, and sphere studies in one
invocation and then creates `dist_lambda_infty.png`.

Run with Rscript from any working directory by supplying the script path.
Outputs default to this directory. For reduced validation, set
`PME_SMOKE_TEST=1` and `PME_OUTPUT_DIR` to a temporary directory. Reduced
settings are not publication results. Full asymptotic runs can take substantial
time. The interval shape-distance calculation also requires `fdasrvf`.
