# Latest-draft alignment audit — October 29, 2026

Manuscript: **Principal Manifold Estimation with Riemannian Templates**.
Source: the user-supplied `Theoretical_Foundations_of_Principal_Manifold_Estimation_with_Non_Euclidean_Templates (3).zip`.
ZIP SHA-256: `20E283300057847166AA5B9FC4BFFD03D00D32CC5B4D64CD132EC929A2FAFBE2`.
GitHub baseline: `7436c8a` on `main`.

## Changes made

- Updated the manuscript title and figure inventory in the repository docs.
- Replaced Figure 5.1 with the exact PNG in the latest ZIP and integrated its
  Unicode headings into the standalone R script. Retained the script's
  output-path setup: the ZIP's copy calls `png(filename = output_file)` without
  defining `output_file`.
- Corrected downstream lambda selection to average the **fold MSEs**:
  `mean(validation_rmse_raw^2)`. Averaging fold RMSEs, or squaring their mean,
  is a different criterion. The inner bandwidth minimizer is unaffected by
  using RMSE versus MSE on the same leave-one-out prediction set.
- Updated all embedded downstream generators consistently. Kept explicitly
  named RMSE diagnostics, and added MSE means and standard errors. Recomputed
  cached summary statistics from the original fold errors; responses,
  predictions, PME functions, and selected bandwidths were retained.
- Re-rendered Figure 5.2 with the correct MSE criterion and axis labels. This
  figure is intentionally different from the still-RMSE PNG in the latest ZIP;
  use the repository's updated PNG in Overleaf to match the draft's definition.
- Restored the production surface-study sample size N=2500, SD=0.02, a
  100-value lambda grid, and the production iteration/tolerance settings in
  both the standalone generator and the plotter's embedded generator.
- Applied strict selection below the elbow (`lambda < cutoff`) and kept the
  full diagnostic curves visible. Updated the residual-average axis text to
  match the latest draft's volume-average interpretation.
- Added the draft's angular-cell coupling TikZ figure under
  `manuscript_native/corrugation_angular_cells.tex` and mapped it in the manifest.

## Verified downstream results

The unprojected outer-fold responses are the validation targets. Squaring each
fold's stored RMSE recovers that fold's MSE exactly; the folds are then averaged.

| Study / method | MSE-selected lambda | Mean outer-fold MSE |
|---|---:|---:|
| Curve / local linear | 1e-8 | 0.003654304 |
| Surface / local linear | 7.19685673001151e-6 | 0.03441698 |
| Curve / kernel diagnostic | 1e-8 | 0.004395057 |
| Surface / kernel diagnostic | 1e-4 | 0.05317250 |

Both local-linear minimizers agree with the draft. No expensive refitting was
needed for the publication Figure 5.2: the original final fits remain valid.
Legacy cached final fits are checked on loading; a changed minimizer requires
regeneration instead of silently pairing a new criterion with an old fit.

## Resolving the parameter conflicts

The latest draft was compared with the actual archived production scripts,
the named data generators, and the supplied figure pixels. A draft number was
not used to overwrite an experimentally supported setting without inspection.

### Surface tuning panels

The archived production `cpme_research_material/simulation_scripts/sphere_job.R`
specifies N=2500 and `noise=0.02`; the smaller driver under
`original_hpc_jobs/sphere_job.R` specifies N=250. The official scripts had
incorrectly inherited the smaller driver's settings.

To distinguish them, the flower data-and-true-surface panel was generated with
the original seed, geometry, camera, and SD=0.02 for N=250, 2500, and 5000.
Comparison with panel (v) of the supplied PNG gave:

| N | Pixel MSE | Black-point-mask intersection / union |
|---:|---:|---:|
| 250 | 0.03308347 | 0.01022785 |
| 2500 | 0.002979428 | 0.8781577 |
| 5000 | 0.04875033 | 0.08506115 |

The N=2500 reconstruction also has 120803 black point pixels versus 121947 in
the original panel, compared with only 13371 for N=250. This corroborates the
production source and the draft's N=2500. SD=0.02 means variance **0.0004**;
the draft's **0.004** does not match the production source or the matching
point pattern. The supplied surface-panel PNGs are retained; full production
surface refits have not been rerun or claimed pixel-identical.

### Additional 1D-in-3D examples

The original `(38) pme examples.Rmd` explicitly uses `noise=0.02` in all four
calls. The associated source functions implement:

- Sinusoidal example before rotation: `(t, sin(2*pi*t), cos(2*pi*t))`.
- Helix before rotation: `(cos(4*pi*t), sin(4*pi*t), 3*t)`.
- Flower in 3D: the stated five-petal radius and `sin(3*pi*t)` height.
- Moon: the adjacent-semicircle construction with its deterministic parameter
  grid and `cos(4*pi*t)` height. Nominal helper argument N=500 emits 1000 points.

The script preserves those original formulas and SD=0.02. The draft instead
states a zero third coordinate for the sinusoidal example, a 2*pi helix, and
noise covariance diagonal 0.04. Those are **manuscript corrections**, not
grounds to replace the source-backed examples with different experiments.
For the existing figures, use `cos(2*pi*t)` as the sinusoidal third coordinate,
`4*pi*t` in the helix's first two coordinates, and noise covariance diagonal
**0.0004** in both interval and circle-template descriptions. Uniform latent
sampling applies to the flower; the moon is constructed by its named helper.

### Other manuscript wording corrections

| Draft statement | Source-backed setting / correction |
|---|---|
| Figure 4.1 uses lambda 1e-9 | The original `(33) paper simulations.Rmd` Monte Carlo block and the official reproduction script use **1e-8**, with 100 repetitions per sample size and SD=0.04. Correct the caption to 1e-8. |
| Planar flower examples have variance 1e-4 | The archived 1D driver and current generator pass SD=0.02, i.e. **4e-4**. This is source-verified; no saved original curve-study objects were available for an independent pixel-level noise check. |
| All cashew examples use bend 1.2*pi | Figure 3.1 does; the original Poisson comparison overrides it with **1.5*pi**, N=1000, SD=0.05. State this example-specific override in the Poisson description. |
| Cashew points are simply bent uniform-sphere draws | The actual helper oversamples 5N sphere points and selects N using k-nearest-neighbor volume weights before adding noise. Include this sampling step for an exact mechanism description. |
| R package is named CPME | The installable R package is **CompactPME**; use `library(CompactPME)`. CPME can remain an informal acronym. |

## Validation and scope

- Before updating, 15 of 16 supplied image assets matched the repository
  byte-for-byte; only Figure 5.1 differed.
- After updating, Figure 5.1 matches the latest ZIP exactly. Figure 5.2 is
  deliberately corrected to MSE. The other fourteen image assets, including
  the original Poisson view/shading and the Ipe drawing, are preserved.
- All R files parse. Regression checks explicitly cover a case where mean
  fold RMSE and mean fold MSE choose different lambdas.
- The combined downstream workflow was tested from an empty input directory
  with reduced settings, and Figure 5.2 was rendered from publication inputs.
- The revised curve-CV and surface-panel workflows were tested with reduced
  settings in separate output directories. Reduced results are not published.
- Full 2500-point surface refits and full 3D/asymptotic studies were not rerun.
  The source-backed manuscript corrections above do not claim new numerical
  experiments or alter the supplied manuscript ZIP.

The immutable draft PNGs cannot substantiate every textual parameter by
themselves. Source-verified findings and the stronger surface-image comparison
are distinguished above so their confidence is not overstated.
