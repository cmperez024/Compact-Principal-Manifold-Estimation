# Five-fold extrinsic regression on a two-dimensional flower

This experiment generalizes the fixed-design one-dimensional flower study to
the CPME project's two-dimensional flower in ambient `R^3`. The response
manifold has an `S^2` template, and the observed external predictor is
two-dimensional.

## Data model

Generate `Z_i = (Z_i1, Z_i2)` independently from `Uniform([-1,1]^2)` and then
condition on the realized design. Define

```text
q_i = [0.5 + 0.45 sin{pi Z_i1 + 0.35 sin(pi Z_i2)}] mod 1
s_i = 0.85 sin{pi Z_i2 + 0.25 sin(pi Z_i1)}
T_i = (sqrt(1-s_i^2) cos(2 pi q_i),
       sqrt(1-s_i^2) sin(2 pi q_i),
       s_i) in S^2
```

The noise-free response is obtained from the project function
`../../shared/datasets.R::flower2d3D_func()` with five petals,
`a = 0.30`, and vertical scale `b = 0.50`. The observed response is

```text
X_i = flower2d3D_func(T_i) + epsilon_i,
epsilon_i ~ N_3(0, 0.04^2 I_3).
```

There is no latent phase error. Once the external covariates have been
generated, the only stochastic error is response-level extrinsic noise.

## Five-fold procedure

For each lambda and outer fold:

1. Fit a two-dimensional CPME surface to only the 320 fold-training responses.
2. Project those training responses onto the fold-specific fitted surface.
3. Select an isotropic two-dimensional regression bandwidth by leave-one-out
   prediction within the fold-training data.
4. Fit either extrinsic kernel regression or extrinsic local-linear regression
   using the pairs `(Z_i, denoised Y_i)`.
5. Predict the 80 excluded responses and project each ambient prediction back
   onto the fold-specific CPME surface.
6. Score predictions against the original raw excluded responses.

The kernel estimator uses Gaussian weights based on Euclidean distance in the
two-dimensional covariate. The local-linear estimator uses the local design
matrix `(1, Z_i1-z_1, Z_i2-z_2)`. Both are multivariate versions of the
estimators in Lin et al. (2017).

The bandwidth is selected inside each training fold, whereas lambda is selected
from the five-fold out-of-fold errors. No excluded response contributes to its
PME surface, bandwidth selection, or regression fit.

Distances for bandwidth selection and final scoring are ambient Euclidean
distances. Intrinsic geodesic distances on the unknown estimated surface are
not required. Projection of regression predictions uses a dense 1,800-point
Fibonacci grid on the fitted `S^2` template.

## Results

Both regression procedures selected `lambda = 1e-4`. This is an interior
minimum of the expanded grid: prediction error increases substantially at
`lambda = 1e-3` and `lambda = 1e-2`.

| Method | Selected lambda | Mean fold bandwidth | Raw OOF RMSE | Oracle OOF RMSE |
|---|---:|---:|---:|---:|
| Extrinsic kernel | 1e-4 | 0.056 | 0.2486 | 0.2369 |
| Extrinsic local-linear | 1e-4 | 0.062 | 0.2004 | 0.1851 |

The extrinsic local-linear estimator reduces raw out-of-fold RMSE by about
19.4% relative to the extrinsic kernel estimator in this realization.

The selected `lambda = 1e-4` fits reached the 15-iteration cap rather than the
strict relative-fidelity tolerance of `1e-4`. Their fidelity histories decrease
monotonically but slowly; the full-data selected fit's last relative change is
approximately `0.00175`. Consequently, the predictive comparison is valid for
the fitted 15-iteration procedure, but the selected surface should not be
described as a fully converged CPME solution.

## Bandwidth sensitivity

The primary experiment selects bandwidth by leave-one-out prediction inside
each fold-training set. A separate sensitivity analysis holds bandwidth fixed
and selects lambda from the outer five-fold errors.

For extrinsic local-linear regression, every fixed bandwidth from `0.04`
through `0.08` selects the smaller `lambda = 7.1969e-6`. Joint minimization over
the tested `(lambda, bandwidth)` grid selects `lambda = 7.1969e-6` and
`h = 0.06`, with raw RMSE `0.1976`. Thus, the local-linear lambda choice is
sensitive to how the nuisance bandwidth is selected.

For extrinsic kernel regression, bandwidths from `0.03` through `0.06` continue
to select `lambda = 1e-4`. Its joint grid minimum is `lambda = 1e-4`, `h = 0.05`,
with raw RMSE `0.2478`.

## Outputs

- `results/fold-performance.csv`: fold-specific bandwidths and errors.
- `results/lambda-performance.csv`: aggregated five-fold performance.
- `results/pme-fit-status.csv`: CPME iteration and convergence information.
- `results/extrinsic-regression-2d-flower-objects.rds`: reusable simulation,
  tuning, fitted-surface, and prediction objects.
- `results/01-kfold-lambda-performance.png`: lambda-selection curves.
- `results/02-data-and-selected-pme-surfaces.png`: true and selected surfaces.
- `results/03-out-of-fold-predictions.png`: honest out-of-fold predictions.
- `results/lambda-selected-by-fixed-bandwidth.csv`: selected lambda for each
  fixed bandwidth and regression method.
- `results/bandwidth-sensitivity-all-combinations.csv`: aggregated errors for
  the complete `(method, lambda, bandwidth)` grid.
- `analyze-bandwidth-sensitivity.R`: reuses the cached PME fits to perform the
  fixed-bandwidth sensitivity analysis without refitting manifolds.

Run `generate-results.R` from anywhere inside the repository. Fold-specific
PME fits are cached under `cache/` so plotting or interrupted reruns do not need
to repeat completed surface fits.
