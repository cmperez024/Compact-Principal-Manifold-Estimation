# Figure 5.2: downstream tuning experiments

This section contains the complete paper-specific experiment code and saved
results for the five-panel downstream extrinsic-regression figure. The CPME
implementation itself is in `CompactPME/` at the repository root.

## Fast reproduction from saved results

From the project root, run:

```powershell
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" "simulation_reproduction\01_main_text\04_downstream_tuning\code_and_results\2d_flower\render-local-linear-composite-figures.R"
```

This recreates `figure_5_2_downstream_tuning.png`. The script uses its own saved
inputs under `figure_5_2_inputs/`, and includes both the 1D and 2D experiment
algorithms to regenerate either input when missing. It does not read the
results in `code_and_results/1d_sd04/` or execute another simulation script.
Required packages are `CompactPME`, `pracma`, and `plot3D`. Supply an absolute
script path to run from any working directory.

The October 29 draft specifies the mean of fold MSEs, so lambda selection now
uses `mean(validation_rmse_raw^2)` and the criterion panels show MSE. This is
not the square of the mean fold RMSE. The cached fold predictions and final
fits are unchanged: the corrected criterion still selects lambda 1e-8 for
the curve and approximately 7.196857e-6 for the surface. Cached objects include
explicit MSE summary columns; the older RMSE columns remain diagnostic.
Legacy cached summaries are upgraded from their fold errors on loading. If
the MSE minimizer differs from a cached final fit, the script requests a new
experiment run instead of displaying an inconsistent fit.

## Full experiment reruns

The one-dimensional SD 0.04 experiment is launched with:

```powershell
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" "simulation_reproduction\01_main_text\04_downstream_tuning\code_and_results\1d_sd04\generate-results.R"
```

The two-dimensional experiment used in the paper is launched with:

```powershell
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" "simulation_reproduction\01_main_text\04_downstream_tuning\code_and_results\2d_flower\generate-results.R" uniform-noninjective 800
```

Both full runs are computationally expensive. The diagnostic renderers also
include their own regeneration code, so these separate generation commands
are optional. `PME_RESULTS_ROOT` redirects experiment inputs/results;
`PME_OUTPUT_DIR` redirects the combined figure. `PME_SMOKE_TEST=1` reduces
sample sizes and tuning grids for validation and must use separate temporary
results/output directories. It does not reproduce publication results.

## Directory roles

- `figure_5_2_downstream_tuning.png`: canonical combined manuscript figure.
- `code_and_results/`: generators, renderers, and canonical saved results in
  short `1d_base`, `1d_sd04`, and `2d_flower` paths.
- `figure_5_2_inputs/`: inputs owned by the combined figure workflow.
- `shared/datasets.R`: historical data definitions retained for reference;
  the simulation scripts include the definitions they use.
