## Two-panel comparison in the same plot3D style used by the project's
## parametric two-dimensional flower figures.

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
notebook_dir <- file.path(
  manuscript_section, "code_and_results", "2d_flower"
)
result_dir <- file.path(notebook_dir, "results-uniform-n800")
object_path <- file.path(
  result_dir, "extrinsic-regression-2d-flower-objects.rds"
)
output_path <- file.path(
  result_dir, "04-true-vs-local-linear-selected-pme.png"
)

local_library <- file.path(repo_root, "packages")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
if (!requireNamespace("pracma", quietly = TRUE)) stop("pracma is required.")
if (!requireNamespace("plot3D", quietly = TRUE)) stop("plot3D is required.")

cart2sph <- pracma::cart2sph
sph2cart <- pracma::sph2cart
ones <- pracma::ones
eye <- pracma::eye
source(file.path(manuscript_section, "shared", "datasets.R"))
source(file.path(repo_root, "CompactPME", "R", "qm.R"))
source(file.path(repo_root, "CompactPME", "R", "spline2d.R"))

objects <- readRDS(object_path)
observed <- objects$data$response_observed
local_linear <- objects$final$local_linear
selected_lambda <- local_linear$choice$lambda

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

true_flower <- function(t) {
  flower2d3D_func(t, r0 = 1, a = 0.30, petals = 5, b = 0.50)
}

true_surface <- generate_surface(true_flower)
estimated_surface <- generate_surface(local_linear$pme$fit)
limits <- range(
  observed,
  true_surface$X, true_surface$Y, true_surface$Z,
  estimated_surface$X, estimated_surface$Y, estimated_surface$Z
)

draw_panel <- function(surface, title) {
  plot3D::scatter3D(
    observed[, 1], observed[, 2], observed[, 3],
    pch = 16, cex = 0.32, col = "black",
    theta = 40, phi = 40,
    main = title,
    bty = "n", ticktype = "detailed",
    xlim = limits, ylim = limits, zlim = limits
  )
  plot3D::surf3D(
    surface$X, surface$Y, surface$Z,
    add = TRUE,
    colvar = surface$Z,
    alpha = 0.60,
    border = "grey",
    bty = "n", colkey = FALSE,
    xlim = limits, ylim = limits, zlim = limits
  )
}

png(output_path, width = 6400, height = 3200, res = 600)
par(mfrow = c(1, 2), mar = c(0.02, 0.1, 3, 0.1))
draw_panel(true_surface, "(i) True parametric flower")
draw_panel(
  estimated_surface,
  sprintf("(ii) Local-linear-selected PME, lambda = %.2e", selected_lambda)
)
dev.off()

cat(normalizePath(output_path, winslash = "/"), "\n")
