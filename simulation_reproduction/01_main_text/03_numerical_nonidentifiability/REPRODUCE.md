# Reproduce Figure 5.1

Status: **complete**.

The Euclidean manuscript's numerical-nonidentifiability construction is
described in Section 5.1 and Appendix F. From the project root, run:

```powershell
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" "simulation_reproduction\01_main_text\03_numerical_nonidentifiability\code\quantile_corrugated_density_curve_grid.R"
```

The implementation uses only base R. It writes
`figure_5_1_numerical_nonidentifiability.png` beside this guide. The output is
a 6000 by 4800 pixel, 400-dpi four-by-five panel containing the nine observed
density alternatives, their nine latent-support manifolds, and the two target
distribution/manifold panels.
