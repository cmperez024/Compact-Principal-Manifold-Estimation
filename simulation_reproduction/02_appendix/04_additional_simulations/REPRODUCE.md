# Reproduce the additional simulations

- `1d3D_sine_spiral.png` and `1d3D_star_moon.png`: run **Section 6 Figures** in `code/one_dimensional_3d_examples.Rmd`. The sine/spiral graphics-device call is commented in the historical source; uncomment its `png(...)` and `dev.off()` lines to save it.
- `dist_lambda_infty.png`: run the circle, interval, and sphere asymptotic sections through the combined `ggsave` call in `code/asymptotic_simulations.Rmd`.

Load CompactPME before running. The asymptotic simulation evaluates several sample sizes and lambda grids and can take substantial time.
