# Compact Principal Manifold Estimation

This repository contains the CompactPME R package and the reproducible figure
workflows associated with the paper described below.

## Associated paper

This repository is the software and simulation-reproduction companion to
["Theoretical Foundations of Principal Manifold Estimation with Non-Euclidean
Templates"](https://arxiv.org/abs/2604.04272) by Kun Meng and Christopher
Perez (2026, arXiv preprint). The `CompactPME/` package contains the principal
manifold estimation implementation, while `simulation_reproduction/` contains
the code and saved inputs used to reproduce the paper's reported simulations
and figures.

## Repository layout

- `CompactPME/`: installable R package.
- `create-CompactPME.Rmd`: authoritative literate source used to create the
  package files.
- `create-CompactPME.html`: rendered, readable version of the package source.
- `simulation_reproduction/`: code, saved inputs, and only the figures used in
  the current manuscript, grouped in manuscript order.

Start with `simulation_reproduction/MANIFEST.csv` to map manuscript figures to
files, and with `simulation_reproduction/CODE_MANIFEST.csv` to locate their
generating code. Each subsection also includes a `REPRODUCE.md` run guide.

The updated independent R workflows in `simulation_reproduction/` are the
official local reproduction versions. Its README lists the canonical entry
points, and the manifests identify each figure's official workflow. Plotters
include the generation code needed for their own missing inputs. Figure 5.2
owns both studies' saved inputs, and the Poisson workflow retains the original
manuscript view and shading through its own reference surface asset.
Reduced validation outputs are kept outside this official tree.

## Installation

```r
remotes::install_github(
  "cmperez024/Compact-Principal-Manifold-Estimation",
  subdir = "CompactPME"
)
```

For local development from the repository root:

```r
devtools::load_all("CompactPME")
```

## Basic use

Given an `N` by `D` data matrix `X`, fit one-dimensional interval-template
estimates over a smoothing-parameter grid with:

```r
library(CompactPME)
fit <- pme_1d_interval(X, 10^(-10:-2), optimize_lambda = TRUE)
```

## Source policy

`create-CompactPME.Rmd` is the package source of record and uses
[`litr`](https://jacobbien.github.io/litr-project/) to generate files under
`CompactPME/`. Edit the source Rmd and re-knit it to update the generated
package; do not edit the generated package files by hand.
The accompanying `create-CompactPME.html` provides a browsable
rendering of that source.

Most functions were developed by Christopher Perez. Kun Meng contributed the
one-dimensional variance-heterogeneity functions `local_var_1d` and `var_het`.
