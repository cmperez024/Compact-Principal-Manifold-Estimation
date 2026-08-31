# Manuscript figure library

This is the canonical, paper-structured figure library for the current
**Euclidean principal-manifold manuscript**.

## Source and preservation policy

- The files here are aligned to the latest supplied Euclidean manuscript project ZIP. They are the publication snapshot, even when a same-named historical workspace file differs.
- Noncanonical, superseded, and unreferenced material is retained locally in
  the repository's ignored archive and waste-bin folders.
- The installable implementation is at `CompactPME/`, while all paper-specific
  experiments live here.
- `MANIFEST.csv` records each curated asset's manuscript role and original ZIP path.
- Each paper subsection has a `code/` directory and a `REPRODUCE.md` guide beside its figures. `CODE_MANIFEST.csv` records the original source of every copied script.

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

The necessary source code is colocated with the figures. Historical Rmd
notebooks, original HPC scripts, and restart caches are outside the publishable
tree. Full simulation scripts regenerate their caches as needed.
