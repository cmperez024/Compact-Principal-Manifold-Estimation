## Paste-ready local-linear criterion + truth + selected-PME figures for the
## current two-dimensional flower experiment and the preceding one-dimensional
## flower experiment.

options(stringsAsFactors = FALSE)

find_repo_root <- function(start = getwd()) {
  here <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(here, "CompactPME", "R"))) return(here)
    parent <- dirname(here)
    if (identical(parent, here)) stop("Could not locate repository root.")
    here <- parent
  }
}

repo_root <- find_repo_root()
setwd(repo_root)
manuscript_section <- file.path(repo_root, "simulation_reproduction",
  "01_main_text", "04_downstream_tuning")
manuscript_experiments <- file.path(manuscript_section, "code_and_results")
local_library <- file.path(repo_root, "packages")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
if (!requireNamespace("pracma", quietly = TRUE)) stop("pracma is required.")
if (!requireNamespace("plot3D", quietly = TRUE)) stop("plot3D is required.")

cart2sph <- pracma::cart2sph
sph2cart <- pracma::sph2cart
bernoulli <- pracma::bernoulli
ones <- pracma::ones
eye <- pracma::eye
source(file.path(manuscript_section, "shared", "datasets.R"))
source(file.path(repo_root, "CompactPME", "R", "qm.R"))
source(file.path(repo_root, "CompactPME", "R", "spline2d.R"))
source(file.path(repo_root, "CompactPME", "R", "krt.R"))
source(file.path(repo_root, "CompactPME", "R", "spline1d_periodic.R"))

draw_criterion <- function(lambda, raw, selected_lambda,
                           selected_raw, main, bandwidth_text) {
  y_limits <- range(raw)
  plot(
    lambda, raw, type = "b", log = "x", pch = 16,
    lwd = 2.4, col = "#D55E00", ylim = y_limits,
    xlab = expression(lambda), ylab = "Five-fold out-of-fold RMSE",
    main = main,
    cex.main = 1.35, cex.lab = 1.28, cex.axis = 1.16
  )
  abline(v = selected_lambda, col = "grey25", lwd = 2.2, lty = 3)
  points(selected_lambda, selected_raw, pch = 21, cex = 1.45,
         bg = "#D55E00", col = "black", lwd = 1.2)
  legend(
    "topleft",
    c("untouched noisy holdout response",
      sprintf("minimizer: lambda = %.2e", selected_lambda)),
    col = c("#D55E00", "grey25"),
    lty = c(1, 3), pch = c(16, NA), lwd = 2.2,
    bty = "n", cex = 1.02
  )
  mtext(bandwidth_text, side = 3, line = 0.25,
        cex = 0.94, col = "grey35")
}

## Two-dimensional flower --------------------------------------------------
notebook_2d <- file.path(
  manuscript_experiments,
  "2d_flower"
)
result_2d <- file.path(notebook_2d, "results-uniform-n800")
objects_2d <- readRDS(file.path(
  result_2d, "extrinsic-regression-2d-flower-objects.rds"
))
profile_2d <- objects_2d$tuning$lambda_performance
profile_2d <- profile_2d[profile_2d$method == "local_linear", ]
profile_2d <- profile_2d[order(profile_2d$lambda), ]
selected_2d <- profile_2d[profile_2d$selected, ]

generate_surface <- function(f, n_grid = 40L) {
  longitude <- seq(0, 2 * pi, length.out = n_grid)
  elevation <- seq(-pi / 2, pi / 2, length.out = n_grid)
  spherical <- expand.grid(longitude, elevation, 1)
  template <- t(apply(spherical, 1L, sph2cart))
  values <- f(template)
  list(
    X = matrix(values[, 1], nrow = n_grid, byrow = TRUE),
    Y = matrix(values[, 2], nrow = n_grid, byrow = TRUE),
    Z = matrix(values[, 3], nrow = n_grid, byrow = TRUE)
  )
}

true_flower_2d <- function(t) {
  flower2d3D_func(t, r0 = 1, a = 0.30, petals = 5, b = 0.50)
}
true_surface_2d <- generate_surface(true_flower_2d)
fit_surface_2d <- generate_surface(objects_2d$final$local_linear$pme$fit)
observed_2d <- objects_2d$data$response_observed
limits_2d <- range(
  observed_2d,
  true_surface_2d$X, true_surface_2d$Y, true_surface_2d$Z,
  fit_surface_2d$X, fit_surface_2d$Y, fit_surface_2d$Z
)

draw_surface_2d <- function(surface, title) {
  plot3D::scatter3D(
    observed_2d[, 1], observed_2d[, 2], observed_2d[, 3],
    pch = 16, cex = 0.45, col = "black",
    theta = 40, phi = 40, main = title,
    bty = "n", ticktype = "detailed",
    xlim = limits_2d, ylim = limits_2d, zlim = limits_2d,
    cex.main = 1.35, cex.lab = 1.24, cex.axis = 1.14
  )
  plot3D::surf3D(
    surface$X, surface$Y, surface$Z, add = TRUE,
    colvar = surface$Z, alpha = 0.60, border = "grey",
    bty = "n", colkey = FALSE,
    xlim = limits_2d, ylim = limits_2d, zlim = limits_2d
  )
}

## One-dimensional flower --------------------------------------------------
notebook_1d <- file.path(
  manuscript_experiments,
  "1d_sd04"
)
result_1d <- file.path(notebook_1d, "results")
drop_1d <- file.path(result_1d, "drop-in")
objects_1d <- readRDS(file.path(
  result_1d, "fixed-design-extrinsic-regression-objects.rds"
))
profile_1d <- read.csv(file.path(
  drop_1d, "lambda-performance-profile.csv"
))
profile_1d <- profile_1d[profile_1d$method == "local_linear", ]
profile_1d <- profile_1d[order(profile_1d$lambda), ]
selected_1d <- profile_1d[profile_1d$selected, ]
observed_1d <- objects_1d$data$response_observed
truth_1d <- objects_1d$display$true_mean_grid
fit_1d <- objects_1d$final$local_linear$pme$fit(
  seq(0, 1, length.out = 2001L)
)
limits_x_1d <- range(observed_1d[, 1], truth_1d[, 1], fit_1d[, 1])
limits_y_1d <- range(observed_1d[, 2], truth_1d[, 2], fit_1d[, 2])

draw_fit_1d <- function(title) {
  plot(
    observed_1d, asp = 1, pch = 16, cex = 0.72,
    col = adjustcolor("grey25", alpha.f = 0.48),
    xlim = limits_x_1d, ylim = limits_y_1d,
    xlab = expression(X[1]), ylab = expression(X[2]),
    main = title,
    cex.main = 1.35, cex.lab = 1.28, cex.axis = 1.16
  )
  lines(truth_1d, col = "black", lwd = 3.2, lty = 2)
  lines(fit_1d, col = "#0072B2", lwd = 3.4)
  legend(
    "topright",
    c("observed response", "true flower", "selected PME manifold"),
    col = c("grey30", "black", "#0072B2"),
    pch = c(16, NA, NA), lty = c(NA, 2, 1),
    lwd = c(NA, 3.2, 3.4), bty = "n", cex = 1.02
  )
}

draw_criterion_1d <- function(title) {
  draw_criterion(
    profile_1d$lambda, profile_1d$mean_raw_rmse,
    selected_1d$lambda, selected_1d$mean_raw_rmse, title,
    sprintf("mean fold-selected bandwidth = %.3f",
            selected_1d$mean_bandwidth)
  )
}

draw_criterion_2d <- function(title) {
  draw_criterion(
    profile_2d$lambda, profile_2d$mean_raw_rmse,
    selected_2d$lambda, selected_2d$mean_raw_rmse, title,
    sprintf("mean fold-selected bandwidth = %.3f", selected_2d$mean_bandwidth)
  )
}

## Combined 2 x 3 publication figure. The 1D fit spans the upper-right cells.
## Write only the figure referenced by the manuscript; component renderings are
## deliberately not emitted into the reproduction tree.
paper_output_combined <- file.path(
  manuscript_section, "figure_5_2_downstream_tuning.png"
)
png(paper_output_combined, width = 5400, height = 3600, res = 300)
layout(
  matrix(c(1, 2, 2,
           3, 4, 5), nrow = 2, byrow = TRUE),
  widths = c(1, 1, 1), heights = c(1, 1)
)
par(mar = c(5.0, 5.2, 4.0, 1.1), mgp = c(3.0, 0.9, 0), cex = 1.08)
draw_criterion_1d("(a) 1D local-linear selection criterion")
par(mar = c(5.0, 5.2, 4.0, 1.1), mgp = c(3.0, 0.9, 0), cex = 1.08)
draw_fit_1d(sprintf(
  "(b) 1D true flower and selected PME, lambda = %.0e",
  selected_1d$lambda
))
par(mar = c(5.0, 5.2, 4.0, 1.1), mgp = c(3.0, 0.9, 0), cex = 1.08)
draw_criterion_2d("(c) 2D local-linear selection criterion")
par(mar = c(0.1, 0.2, 4.0, 0.2), cex = 1.08)
draw_surface_2d(true_surface_2d, "(d) 2D true parametric flower")
draw_surface_2d(
  fit_surface_2d,
  sprintf("(e) 2D selected PME, lambda = %.2e", selected_2d$lambda)
)
dev.off()

cat(normalizePath(paper_output_combined, winslash = "/"), "\n")
