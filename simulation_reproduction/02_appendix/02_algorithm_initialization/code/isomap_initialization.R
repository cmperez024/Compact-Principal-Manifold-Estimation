# Reproduce iso_1d_2d.png from raw simulated data.

if (!requireNamespace("CompactPME", quietly = TRUE)) {
  stop("Install the CompactPME package before running this script.")
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L])))
} else {
  normalizePath(getwd())
}
output_file <- Sys.getenv(
  "PME_ISOMAP_OUTPUT",
  unset = file.path(script_dir, "..", "iso_1d_2d.png")
)

moon_curve <- function(n, inner_radius = 1, outer_radius = 2, noise_sd = 0) {
  stopifnot(outer_radius > inner_radius)
  n_main <- round(4 * n / 5)
  n_cap <- round(n / 5)
  cap_radius <- (outer_radius - inner_radius) / 2
  cap_center <- (inner_radius + outer_radius) / 2
  theta_outer <- seq(pi / 2, 3 * pi / 2, length.out = n_main)
  phi_bottom <- seq(-pi / 2, pi / 2, length.out = n_cap)
  theta_inner <- seq(3 * pi / 2, pi / 2, length.out = n_main)
  phi_top <- seq(-pi / 2, pi / 2, length.out = n_cap)
  signal <- cbind(
    c(outer_radius * cos(theta_outer), cap_radius * cos(phi_bottom),
      inner_radius * cos(theta_inner), cap_radius * cos(phi_top)),
    c(outer_radius * sin(theta_outer),
      -cap_center + cap_radius * sin(phi_bottom),
      inner_radius * sin(theta_inner),
      cap_center + cap_radius * sin(phi_top))
  )
  signal + matrix(rnorm(length(signal), sd = noise_sd), ncol = 2L)
}

set.seed(1)
X <- moon_curve(300L, noise_sd = 0.02)
true_curve <- moon_curve(1000L)
evaluation_grid <- seq(0, 1, length.out = 1001L)

isomap_1d <- vegan::isomap(stats::dist(X), ndim = 1L, k = 10L)
fit_1d <- CompactPME::spline1d_periodic(
  X, as.numeric(isomap_1d$points), lambda = 1e-6
)

isomap_2d <- vegan::isomap(stats::dist(X), ndim = 2L, k = 10L)
circle_index <- CompactPME:::project_circle(isomap_2d$points)
fit_2d <- CompactPME::spline1d_periodic(X, circle_index, lambda = 1e-6)

grDevices::png(output_file, width = 8, height = 4, units = "in", res = 600)
graphics::par(mfrow = c(1, 2), mar = c(3.2, 3.2, 1.2, 0.7))
true_color <- grDevices::adjustcolor("grey20", alpha.f = 0.5)

graphics::plot(true_curve, type = "l", asp = 1, col = true_color, lwd = 2,
               xlab = "x", ylab = "y")
graphics::lines(fit_1d(evaluation_grid), col = "darkgreen", lwd = 2)

graphics::plot(true_curve, type = "l", asp = 1, col = true_color, lwd = 2,
               xlab = "x", ylab = "y")
graphics::lines(fit_2d(evaluation_grid), col = "blue", lty = 2, lwd = 2)
grDevices::dev.off()

message("Wrote ", output_file)
