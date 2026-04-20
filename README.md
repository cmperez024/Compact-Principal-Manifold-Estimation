# Compact Principal Manifold Estimation

## Installation and Getting Started

To install the package, please use the commands below. Note that you must specify the subdir for it to work correctly.

``` r
library(devtools)
install_github("cmperez024/Compact-Principal-Manifold-Estimation", subdir="CompactPME")
```

To use the package, you must first have an $N \times D$ data matrix where $D$ is the ambient dimension, and some smoothing parameter value $\lambda$. The code

``` r
# Given some data, X
fit <- pme_1d_interval(X, 10^(-10:-2), optimize_lambda=T)
```

will return a manifold estimates for each of the $\lambda$ values in $\{10^{-10}, \ldots, 10^{-2}\}.$ It will also perform the optimization routine, where one may observe the coefficient of heterogeneity values for each fit.

## Files

Outside of the directory CompactPME (which contains the project files), we also have the source code for the package which is written using the [litr](https://jacobbien.github.io/litr-project/) package in R. This is available in its raw Rmd form as create-CompactPME.Rmd. In the Simulation Scripts folder, there are scripts that can be run independently of the package to reproduce results in the paper. Also in this folder is the datasets file which contains multiple functions for generating simulated data.

## Contributions

Most functions developed by Christopher Perez. Code for variance heterogeneity in 1d (local_var_1d and var_het) provided by Kun Meng.

## Citing this Package

This package is associated with the paper ["Theoretical Foundations of Principal Manifold Estimation with Non-Euclidean Templates" (Meng, Perez, 2026, arXiv Pre-Print)](https://arxiv.org/abs/2604.04272).
