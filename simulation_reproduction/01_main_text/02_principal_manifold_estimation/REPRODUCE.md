# Reproduce the principal-manifold figures

- `2d3D_star_moon.png`: run the **Section 3 figures** chunks in `code/surface_and_curve_examples.Rmd`.
- `PME_monte_carlo_star2.png`: run **Demo Monte Carlo** in `code/paper_simulations_and_monte_carlo.Rmd`. It performs 100 fits for each of `N = 100, 250, 500, 1000`.

Load the preserved CompactPME implementation before running the notebook sections:

```r
devtools::load_all("CompactPME")
```

These are computationally expensive. Write reruns to candidate filenames and compare them before replacing the publication PNGs.
