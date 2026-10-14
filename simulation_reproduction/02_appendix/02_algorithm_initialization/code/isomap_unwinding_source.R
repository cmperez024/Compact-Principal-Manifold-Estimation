# Reproduce the ISOMAP unwinding diagnostic and the initial-versus-converged
# PME comparison from the historical simulation notebook.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1L])))
figure_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
repo_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."),
                           mustWork = TRUE)
local_library <- file.path(repo_root, "packages")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))

required_packages <- c("CompactPME", "vegan", "scico")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install the required packages: ", paste(missing_packages, collapse = ", "))
}

unwind_output <- Sys.getenv(
  "PME_ISOMAP_UNWIND_OUTPUT",
  unset = file.path(figure_dir, "isomap_unwind.png")
)
comparison_output <- Sys.getenv(
  "PME_ISOMAP_COMPARISON_OUTPUT",
  unset = file.path(figure_dir, "isomap_unwind_first_vs_converged.png")
)

moon_curve <- function(n, inner_radius = 1, outer_radius = 2, noise_sd = 0) {
  stopifnot(outer_radius > inner_radius)
  n_main <- round(4 * n / 5)
  n_cap <- round(n / 5)
  cap_radius <- (outer_radius - inner_radius) / 2
  cap_center <- (inner_radius + outer_radius) / 2
  theta_outer <- seq(pi / 2, 3 * pi / 2, length.out = n_main)
  theta_inner <- seq(3 * pi / 2, pi / 2, length.out = n_main)
  phi <- seq(-pi / 2, pi / 2, length.out = n_cap)
  signal <- cbind(
    c(outer_radius * cos(theta_outer), cap_radius * cos(phi),
      inner_radius * cos(theta_inner), cap_radius * cos(phi)),
    c(outer_radius * sin(theta_outer),
      -cap_center + cap_radius * sin(phi),
      inner_radius * sin(theta_inner),
      cap_center + cap_radius * sin(phi))
  )
  signal + matrix(stats::rnorm(length(signal), sd = noise_sd), ncol = 2L)
}

set.seed(1)
X <- moon_curve(300L, noise_sd = 0.06)
embedding <- vegan::isomap(stats::dist(X), ndim = 2L, k = 10L)$points
# Eigenvector signs are arbitrary; fix the horizontal orientation to the
# orientation used in the manuscript snapshot.
embedding[, 1L] <- -embedding[, 1L]
point_colors <- scico::scico(nrow(X), palette = "romaO")

grDevices::png(unwind_output, width = 875, height = 540)
graphics::par(mfrow = c(1, 2))
graphics::plot(X, col = point_colors, xlab = "x", ylab = "y")
graphics::plot(embedding, col = point_colors,
               xlab = "Dim1", ylab = "Dim2")
grDevices::dev.off()

fit_curve <- function(use_isomap, iterations) {
  set.seed(1)
  CompactPME::pme_1d_periodic(
    X,
    lambdas = 1e-8,
    init_iso = use_isomap,
    k_init = 10L,
    max_iter = iterations,
    optimize_lambda = FALSE,
    parallel = FALSE,
    approx = FALSE
  )$spline_list[[1L]]
}

naive_initial <- fit_curve(FALSE, 1L)
isomap_initial <- fit_curve(TRUE, 1L)
naive_converged <- fit_curve(FALSE, 10L)
isomap_converged <- fit_curve(TRUE, 10L)

evaluation_grid <- seq(0, 1, length.out = 1000L)
curves <- list(
  naive_initial(evaluation_grid),
  isomap_initial(evaluation_grid),
  naive_converged(evaluation_grid),
  isomap_converged(evaluation_grid)
)
x_limits <- range(vapply(curves, function(x) range(x[, 1L]), numeric(2L)))
y_limits <- range(vapply(curves, function(x) range(x[, 2L]), numeric(2L)))

grDevices::png(comparison_output, width = 8, height = 4,
               units = "in", res = 600)
graphics::par(mfrow = c(1, 2), mar = c(3.2, 3.2, 1.2, 0.7))
graphics::plot(curves[[1L]], type = "l", asp = 1, col = "red", lwd = 2,
               xlim = x_limits, ylim = y_limits, xlab = "x", ylab = "y")
graphics::lines(curves[[2L]], col = "blue", lty = 2, lwd = 2)
graphics::plot(curves[[3L]], type = "l", asp = 1, col = "red", lwd = 2,
               xlim = x_limits, ylim = y_limits, xlab = "x", ylab = "y")
graphics::lines(curves[[4L]], col = "blue", lty = 2, lwd = 2)
grDevices::dev.off()

message("Wrote ", normalizePath(unwind_output, winslash = "/"))
message("Wrote ", normalizePath(comparison_output, winslash = "/"))
