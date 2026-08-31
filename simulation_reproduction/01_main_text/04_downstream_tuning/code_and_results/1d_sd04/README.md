# Fixed-design extrinsic regression with response-noise SD 0.04

This repeats the `../1d_base/` experiment with

\[
\varepsilon_i\overset{\mathrm{iid}}{\sim}N_2(0,0.04^2I_2)
\]

instead of variance `0.08^2`. The seed, fixed covariates, deterministic phase
map, standardized Gaussian error realization, five folds, PME lambda grid,
and regression-bandwidth grid are unchanged. Consequently, differences from
the SD 0.08 experiment isolate the effect of halving response-noise SD.

Run `generate-results.R` from the repository root. Use `plot-selected-pme.R`
to recreate the final figure from the saved RDS without refitting.

## Result

All 55 fold-specific PME fits converged. For each candidate `lambda` and outer
fold, the bandwidth was selected by leave-one-out prediction of the projected
training responses; untouched outer-fold responses were used only to select
`lambda`. Both downstream criteria selected `lambda = 1e-8`. For the kernel
estimator, the mean fold-selected bandwidth was `0.0064`, the full-data
bandwidth was `0.006`, and the raw and simulation-oracle RMSEs were `0.06626`
and `0.03483`. For the local-linear estimator, the corresponding values were
`0.0092`, `0.008`, `0.06044`, and `0.02215`.
