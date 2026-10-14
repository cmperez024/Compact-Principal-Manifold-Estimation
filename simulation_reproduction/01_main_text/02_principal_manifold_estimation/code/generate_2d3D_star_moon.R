# Reproduce the 2 x 4 star/moon PME surface figure and cache its fits.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1L])))
figure_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
repo_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."), mustWork = TRUE)
local_library <- file.path(repo_root, "packages")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))

required_packages <- c("CompactPME", "FNN", "plot3D", "pracma")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the required packages before running: ",
    paste(missing_packages, collapse = ", ")
  )
}

lambdas <- 10^c(-12, -5.5, -3, -2.2)
fit_path <- Sys.getenv(
  "PME_STAR_MOON_RDS",
  file.path(figure_dir, "2d3D_star_moon_fits.rds")
)
figure_path <- Sys.getenv(
  "PME_STAR_MOON_OUTPUT",
  file.path(figure_dir, "2d3D_star_moon.png")
)
refit <- tolower(Sys.getenv("PME_STAR_MOON_REFIT", "true")) %in%
  c("1", "true", "yes")

# The historical fitter lists Riemann as a foreach worker package although the
# present optimize_lambda = FALSE run never uses it. Keep this compatibility
# adjustment local to this script instead of adding an unrelated dependency.
fit_surface <- CompactPME::pme_2d_periodic
fit_surface_body <- paste(deparse(body(fit_surface)), collapse = "\n")
fit_surface_body <- gsub(
  'c("stats", "pracma", "Riemann")',
  'c("stats", "pracma")',
  fit_surface_body,
  fixed = TRUE
)
body(fit_surface) <- parse(text = fit_surface_body)[[1L]]
environment(fit_surface) <- asNamespace("CompactPME")

random_sphere <- function(n) {
  u <- runif(n, -1, 1)
  azimuth <- runif(n, 0, 2 * pi)
  polar <- acos(u)
  cbind(
    sin(polar) * cos(azimuth),
    sin(polar) * sin(azimuth),
    cos(polar)
  )
}

flower_2d3d <- function(n = 300L, petals = 6L, noise = 0.01) {
  theta <- runif(n, 0, 2 * pi)
  z_unit <- runif(n, -1, 1)
  radius_xy <- 1 + 0.3 * cos(petals * theta)
  signal <- cbind(
    radius_xy * sqrt(1 - z_unit^2) * cos(theta),
    radius_xy * sqrt(1 - z_unit^2) * sin(theta),
    0.5 * z_unit
  )
  signal + matrix(rnorm(3L * n, sd = noise), nrow = n, ncol = 3L)
}

moon_2d3d <- function(n = 300L, oversample = 5L, bend = 1.2 * pi,
                      radius = 1, bend_radius = 2, noise = 0.01) {
  sphere <- random_sphere(n * oversample)
  x <- radius * sphere[, 1L]
  y <- radius * sphere[, 2L]
  z <- radius * sphere[, 3L]
  points <- cbind(
    cos(bend * x / (2 * radius)) * (bend_radius * radius + y),
    sin(bend * x / (2 * radius)) * (bend_radius * radius + y),
    z
  )
  kth_distance <- FNN::get.knn(points, k = 10L)$nn.dist[, 10L]
  keep <- sample(seq_len(nrow(points)), n, prob = kth_distance^3)
  points[keep, , drop = FALSE] +
    matrix(rnorm(3L * n, sd = noise), nrow = n, ncol = 3L)
}

if (refit || !file.exists(fit_path)) {
  set.seed(1)
  star_data <- flower_2d3d()
  star_fit <- fit_surface(
    star_data,
    lambdas = lambdas,
    parallel = FALSE,
    optimize_lambda = FALSE
  )
  moon_data <- moon_2d3d()
  moon_fit <- fit_surface(
    moon_data,
    lambdas = lambdas,
    parallel = FALSE,
    optimize_lambda = FALSE,
    init_iso = TRUE
  )
  fit_bundle <- list(
    schema_version = 1L,
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    seed = 1L,
    lambdas = lambdas,
    settings = list(
      n = 300L,
      noise = 0.01,
      star_petals = 6L,
      moon_bend_radius = 2,
      moon_bend = 1.2 * pi,
      moon_init_iso = TRUE
    ),
    data = list(star = star_data, moon = moon_data),
    fits = list(star = star_fit, moon = moon_fit)
  )
  saveRDS(fit_bundle, fit_path, compress = "xz")
} else {
  fit_bundle <- readRDS(fit_path)
}

surface_grid <- function(fit, grid_size = 50L) {
  longitude <- seq(0, 2 * pi, length.out = grid_size)
  latitude <- seq(-pi / 2, pi / 2, length.out = grid_size)
  spherical_grid <- expand.grid(longitude, latitude, 1)
  cartesian_grid <- t(apply(spherical_grid, 1L, pracma::sph2cart))
  values <- fit(cartesian_grid)
  list(
    x = matrix(values[, 1L], grid_size, grid_size, byrow = TRUE),
    y = matrix(values[, 2L], grid_size, grid_size, byrow = TRUE),
    z = matrix(values[, 3L], grid_size, grid_size, byrow = TRUE)
  )
}

plot_panel <- function(data, fit, theta, phi, lambda) {
  surface <- surface_grid(fit)
  limits <- range(data)
  plot3D::scatter3D(
    data[, 1L], data[, 2L], data[, 3L],
    pch = 16, cex = 0.42, col = "black",
    theta = theta, phi = phi,
    main = "",
    xlab = "x", ylab = "y", zlab = "z",
    bty = "g", ticktype = "detailed",
    xlim = limits, ylim = limits, zlim = limits
  )
  plot3D::surf3D(
    surface$x, surface$y, surface$z,
    add = TRUE, colvar = surface$z,
    alpha = 0.6, border = "grey", colkey = FALSE,
    bty = "g", xlim = limits, ylim = limits, zlim = limits
  )
  lambda_label <- bquote(lambda == .(formatC(lambda, format = "e", digits = 2)))
  legend(
    "top",
    legend = lambda_label,
    inset = c(0, 0.065),
    bty = "o",
    bg = "white",
    box.col = "white",
    cex = 0.98,
    text.font = 2,
    x.intersp = 0.3,
    y.intersp = 0.95
  )
}

png(figure_path, width = 7200, height = 3300, res = 600)
on.exit(dev.off(), add = TRUE)
par(
  mfrow = c(2, 4),
  mar = c(0.25, 0.35, 0.25, 0.20),
  oma = c(0, 0, 0, 0),
  mgp = c(0.55, 0.12, 0),
  tcl = -0.16,
  xpd = FALSE
)
for (i in seq_along(lambdas)) {
  plot_panel(
    fit_bundle$data$star,
    fit_bundle$fits$star$spline_list[[i]],
    theta = 20, phi = 50,
    lambda = lambdas[i]
  )
}
for (i in seq_along(lambdas)) {
  plot_panel(
    fit_bundle$data$moon,
    fit_bundle$fits$moon$spline_list[[i]],
    theta = 200, phi = 40,
    lambda = lambdas[i]
  )
}
dev.off()
on.exit(NULL, add = FALSE)

message("Saved fits: ", normalizePath(fit_path, mustWork = TRUE))
message("Saved figure: ", normalizePath(figure_path, mustWork = TRUE))
