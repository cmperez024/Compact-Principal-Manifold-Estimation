# Figure 5.2: downstream tuning experiments

This section contains the complete paper-specific experiment code and saved
results for the five-panel downstream extrinsic-regression figure. The CPME
implementation itself is in `CompactPME/` at the repository root.

## Fast reproduction from saved results

From the project root, run:

```powershell
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" "simulation_reproduction\01_main_text\04_downstream_tuning\code_and_results\2d_flower\render-local-linear-composite-figures.R"
```

This recreates `figure_5_2_downstream_tuning.png`. Intermediate one- and
two-dimensional renderings are not written, so the reproduction tree continues
to contain only manuscript figures. The command uses saved RDS and CSV objects,
so it does not refit the manifolds. It requires the R packages `pracma` and
`plot3D`.

## Full experiment reruns

The one-dimensional SD 0.04 experiment is launched with:

```powershell
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" "simulation_reproduction\01_main_text\04_downstream_tuning\code_and_results\1d_sd04\generate-results.R"
```

The two-dimensional experiment used in the paper is launched with:

```powershell
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" "simulation_reproduction\01_main_text\04_downstream_tuning\code_and_results\2d_flower\generate-results.R" uniform-noninjective 800
```

Both full runs are computationally expensive. They create restart caches beside
the scripts; those generated caches are intentionally excluded from Git. The
canonical final saved objects needed by the fast renderer remain included.

## Directory roles

- `figure_5_2_downstream_tuning.png`: canonical combined manuscript figure.
- `code_and_results/`: generators, renderers, and canonical saved results in
  short `1d_base`, `1d_sd04`, and `2d_flower` paths.
- `shared/datasets.R`: simulation data functions required by the 2D scripts.
