# Reproduce the principal-manifold figures

- `2d3D_star_moon.png`: run the **Section 3 figures** chunks in
  `code/surface_and_curve_examples.Rmd`.
- `PME_monte_carlo_star2.png`: from the repository root, run:

```powershell
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" "simulation_reproduction\01_main_text\02_principal_manifold_estimation\code\generate_figure_4_1_consistency_study.R"
```

The script contains only the Figure 4.1 experiment and renderer. It performs
100 fits for each of `N = 100, 250, 500, 1000`, using noise standard deviation
`0.04` and fixed `lambda = 1e-8`. Set `PME_FIGURE_OUTPUT` to write a candidate
figure elsewhere instead of replacing the canonical PNG.

Install CompactPME first. For interactive work, load the local implementation
with:

```r
devtools::load_all("CompactPME")
```

The consistency study is computationally expensive because it performs 400 PME
fits. Compare a candidate rerun with the canonical figure before replacing it.
