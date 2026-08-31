# Compact Principal Manifold Estimation

This repository contains the CompactPME R package and the reproducible figure
workflows for its associated Euclidean principal-manifold manuscript.

## Repository layout

- `CompactPME/`: installable R package.
- `create-CompactPME.Rmd`: authoritative literate source used to create the
  package files.
- `simulation_reproduction/`: code, saved inputs, and only the figures used in
  the current manuscript, grouped in manuscript order.

Start with `simulation_reproduction/MANIFEST.csv` to map manuscript figures to
files, and with `simulation_reproduction/CODE_MANIFEST.csv` to locate their
generating code. Each subsection also includes a `REPRODUCE.md` run guide.

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
fit <- pme_1d_interval(X, 10^(-10:-2), optimize_lambda = TRUE)
```

## Source policy

`create-CompactPME.Rmd` is the package source of record and uses
[`litr`](https://jacobbien.github.io/litr-project/) to generate files under
`CompactPME/`. Changes to generated package code should also be reflected in
the source Rmd.

Most functions were developed by Christopher Perez. Kun Meng contributed the
one-dimensional variance-heterogeneity functions `local_var_1d` and `var_het`.
