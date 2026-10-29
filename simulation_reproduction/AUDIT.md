# Simulation reproduction audit

Audit target: `simulation_reproduction/`
Independence pass completed: October 13, 2026.

Latest-draft follow-up: see `DRAFT_ALIGNMENT_2026-10-29.md`. That audit
supersedes the historical Figure 5.2 RMSE rendering and the older surface
driver settings. Figure 5.2 now uses mean fold MSE; the verified production
surface configuration is N=2500 and SD=0.02.

The independent workflows are now the official local versions. The repository
README and simulation README identify them as the maintained source of record;
the manifests map each publication asset to its canonical entry point.
Optional generators remain available for precomputation. Reduced baseline
test results and the superseded Figure 3.1 duplicate filename were moved
recoverably into the ignored validation archive before promotion.

## Reproduction policy

Canonical simulations are plain R scripts. Each figure workflow includes the
data generation and calculation code needed to regenerate its own missing
inputs. It obtains PME fitting and kernel functions from the installed
CompactPME package; it does not source another simulation or package R file.
Other R packages may be imported. Run scripts with Rscript using their paths;
paths are resolved from the script rather than the caller's working directory.

Saved inputs can accelerate a run but are optional. Figure 5.2 owns its 1D and
2D input objects under `figure_5_2_inputs/` and includes both experiment
algorithms. Its renderer no longer reads `code_and_results/1d_sd04/` results.
Separate generators remain optional tools for precomputing the same study.

The star-convex illustration is an Ipe drawing. Its supplied manuscript PNG
is retained exactly and it has no R generator.

## Changes in this pass

- Resolved all twelve downstream scripts' repository paths from their files.
- Added embedded regeneration to the overview, CV, sphere-panel, downstream
  diagnostic, and combined-figure plotters.
- Removed copied package fitting/kernel implementations from the older 1D
  and sphere generators, replacing them with installed package bindings.
- Enabled the flower TDA calculation required by the sphere-panel workflow.
- Extracted the two 1D-in-3D figures from the historical multi-experiment code,
  embedding the necessary data definitions and using current CompactPME APIs.
- Extracted the three large-lambda studies into one independent asymptotic
  workflow that computes its own inputs rather than reading other experiments.
- Automated the Poisson calculation and panel combination. After the view
  correction, the default uses this workflow's original shaded snapshot asset
  to preserve the unrecorded manual camera exactly. Fresh rgl rendering is
  available explicitly; the flat orthographic renderer has been removed.
- Replaced the stale canonical Figure 3.1 PNG with the manuscript-matching
  7200 x 3300 rendering.
- Removed the unused Riemann worker dependency in create-CompactPME.Rmd,
  regenerated the package with litr, and installed the generated package.

## Validation

All 27 R files parse, no R Markdown files remain in the reproduction folder,
and no source() calls remain. Both manifests resolve to existing local files.

Executed from a working directory outside the project:

| Workflow | Validation |
|---|---|
| Figure 5.2 combined renderer | Publication cached inputs: output byte-for-byte identical to manuscript PNG. Fresh empty input directory: reduced experiment generated both studies and the final figure. |
| Overview plotter | Fresh directory: reduced generation, persistence calculation, PME fits, and all output figures completed. |
| Sphere panels | Reduced generation of both persistence diagrams and PME fits completed; both final panels assembled successfully. |
| 1D-in-3D appendix figures | Both figures generated successfully with reduced sample/iteration settings. |
| Large-lambda appendix figure | All three studies and combined figure completed with reduced settings. |
| SD 0.04 diagnostic plotter | Fresh directory: reduced generation and diagnostics completed. |
| Baseline 1D plotter | Fresh directory: reduced generation and final comparison completed. |
| Fixed-bandwidth 2D plotter | Fresh directory: reduced fitting, sensitivity calculation, and all three static figures completed. |
| Interactive fixed-bandwidth plotter | Fresh directory: reduced fitting and HTML output completed. |
| Poisson comparison | Full 1000-observation calculation completed. Corrected default combined output matches the original manuscript image pixel-for-pixel. Fresh shaded rgl rendering was also tested. |

Reduced tests use PME_SMOKE_TEST=1 and separate local validation directories.
They verify execution and file independence, not publication numerical results
or full convergence. Full expensive overview, sphere, downstream, 1D-in-3D,
and asymptotic reruns were not executed. These tests must not overwrite or
replace publication inputs.

## Appearance and dependencies

At the October 13 audit, Figure 5.2 and the corrected Figure 3.1 matched the
then-supplied manuscript assets. The October 29 follow-up re-renders Figure 5.2
in MSE units and updates Figure 5.1 from the new draft. The original Ipe
star-convex PNG is unchanged. The default
Poisson figure now preserves the historical view and shading exactly through
its own original surface snapshot. The optional fresh rgl mode does not claim
an exact camera or triangulation match to that historical snapshot.
The extracted appendix workflows now have explicit fixed seeds; full
publication-size visual comparisons remain to be performed for those figures.

SurfaceReconstruction built successfully using C++17 with current RcppCGAL.
Its installation instructions are recorded in that section's REPRODUCE.md.
The previous installed CompactPME package is preserved in
packages/CompactPME.pre-independence-backup.

This report records local validation. Publication of the official workflows
is recorded in the repository's Git history; reduced validation outputs
remain local and ignored.
