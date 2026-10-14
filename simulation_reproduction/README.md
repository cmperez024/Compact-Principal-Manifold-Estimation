# Manuscript figure library

This is the canonical, paper-structured figure library for the current
**Euclidean principal-manifold manuscript**.

The updated independent R scripts in this directory are the official local
reproduction versions, effective October 13, 2026. Use the entry points below.
The original notebooks and archive files record provenance; the R scripts
listed here are the maintained workflows.

## Official entry points

| Figure or study | Official R entry point |
|---|---|
| Figure 2.1 overview | [plot_1d_figures.R](01_main_text/01_overview_and_preliminaries/code/plot_1d_figures.R) |
| Figure 3.1 star/moon surfaces | [generate_2d3D_star_moon.R](01_main_text/02_principal_manifold_estimation/code/generate_2d3D_star_moon.R) |
| Figure 4.1 consistency study | [generate_figure_4_1_consistency_study.R](01_main_text/02_principal_manifold_estimation/code/generate_figure_4_1_consistency_study.R) |
| Figure 5.1 nonidentifiability | [quantile_corrugated_density_curve_grid.R](01_main_text/03_numerical_nonidentifiability/code/quantile_corrugated_density_curve_grid.R) |
| Figure 5.2 downstream tuning | [render-local-linear-composite-figures.R](01_main_text/04_downstream_tuning/code_and_results/2d_flower/render-local-linear-composite-figures.R) |
| Appendix Poisson comparison | [poisson_comparison.R](02_appendix/01_surface_reconstruction_comparison/code/poisson_comparison.R) |
| Appendix 1D/2D ISOMAP | [isomap_initialization.R](02_appendix/02_algorithm_initialization/code/isomap_initialization.R) |
| Appendix ISOMAP unwinding | [isomap_unwinding_source.R](02_appendix/02_algorithm_initialization/code/isomap_unwinding_source.R) |
| Appendix 1D tuning | [plot_1d_tuning_figures.R](02_appendix/03_tuning_parameter_selection/code/plot_1d_tuning_figures.R) |
| Appendix sphere tuning panels | [plot_sphere_panels.R](02_appendix/03_tuning_parameter_selection/code/plot_sphere_panels.R) |
| Appendix 1D curves in 3D | [one_dimensional_3d_examples.R](02_appendix/04_additional_simulations/code/one_dimensional_3d_examples.R) |
| Appendix large-lambda limits | [asymptotic_simulations.R](02_appendix/04_additional_simulations/code/asymptotic_simulations.R) |

Run these with Rscript using the script path from any working directory.
Each workflow generates its own missing simulation inputs and imports the
installed CompactPME package. The separate generation scripts remain official
optional precomputation tools for their respective studies.

The Ipe star-convex illustration and manuscript-native TikZ figures retain
their original drawing sources. For the Poisson publication figure, include
`02_appendix/01_surface_reconstruction_comparison/assets/poisson_surface_manuscript.png`:
it preserves the original manual view and shading. Fresh rgl rendering is an
explicit alternative described in that section's run guide.

## Source and preservation policy

- The files here are aligned to the latest supplied Euclidean manuscript project ZIP. They are the publication snapshot, even when a same-named historical workspace file differs.
- Noncanonical, superseded, and unreferenced material is retained locally in
  the repository's ignored archive and waste-bin folders.
- The installable implementation is at `CompactPME/`, while all paper-specific
  experiments live here.
- `MANIFEST.csv` records each curated asset's manuscript role and original ZIP path.
- Each paper subsection has a `code/` directory and a `REPRODUCE.md` guide beside its figures. `CODE_MANIFEST.csv` records provenance and identifies official entry points and optional precomputation tools. `MANIFEST.csv` maps assets to their official entry points.

## Folder structure

```text
simulation_reproduction/
|-- 01_main_text/
|   |-- 01_overview_and_preliminaries/
|   |-- 02_principal_manifold_estimation/
|   |-- 03_numerical_nonidentifiability/
|   `-- 04_downstream_tuning/
`-- 02_appendix/
|   |-- 01_surface_reconstruction_comparison/
|   |-- 02_algorithm_initialization/
|   |-- 03_tuning_parameter_selection/
    `-- 04_additional_simulations/
```

## Manuscript-native figures

Three figures are drawn directly in LaTeX/TikZ and therefore do not have standalone image files in this library:

1. The fitted curve, projection, and residual illustration (`fig:proj-curve`).
2. The ellipse and medial axis illustration (`fig:ellipse-medial-axis`).
3. The projection-adaptation basin-of-attraction illustration (`fig:placeholder`).

Their source remains in the main manuscript's LaTeX file. They are recorded in
`MANIFEST.csv` as manuscript-native figures and do not require R simulation
code.

## Status note on the current PDF

The current manuscript contains main-text Figures 1.1, 2.1--2.2, 3.1--3.2,
4.1, and 5.1--5.2. Figure 5.1 and its base-R generator are under
`01_main_text/03_numerical_nonidentifiability/`. The five-panel downstream-
tuning Figure 5.2 is under `01_main_text/04_downstream_tuning/` together with
its generators, renderer, and the saved results required for fast rendering.

## Working convention

Use this directory first when reproducing, editing, compiling, or submitting manuscript figures. If a figure is regenerated, compare it to the curated version before replacing the publication snapshot.

## Reproduction status

The maintained simulations are R scripts. Historical Rmd notebooks remain in
the local archive. The retained HPC launcher is a historical submission example.
Full workflows regenerate their own caches as needed, while Figure 5.2's
publication inputs are included under `figure_5_2_inputs/`. Validation outputs
and superseded duplicate filenames are archived under the ignored
`.independence_validation/` directory at the project root. See `AUDIT.md` for
the completed execution checks and outstanding full-size comparisons.
