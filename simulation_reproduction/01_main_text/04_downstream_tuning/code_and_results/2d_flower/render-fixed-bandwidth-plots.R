## Static figures for h = 0.06 and the local-linear selected lambda.

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
result_dir <- file.path(notebook_dir, "results")
local_library <- file.path(repo_root, "packages")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))

cart2sph <- pracma::cart2sph
sph2cart <- pracma::sph2cart
source(file.path(manuscript_section, "shared", "datasets.R"))
source(file.path(repo_root, "CompactPME", "R", "qm.R"))

objects <- readRDS(file.path(
  result_dir, "fixed-bandwidth-local-linear-objects.rds"
))
sensitivity <- read.csv(file.path(
  result_dir, "bandwidth-sensitivity-all-combinations.csv"
))

camera_rotation <- function(azimuth = 35, elevation = 23) {
  az <- azimuth * pi / 180
  el <- elevation * pi / 180
  rz <- matrix(c(cos(az), -sin(az), 0,
                 sin(az),  cos(az), 0,
                 0, 0, 1), 3, 3, byrow = TRUE)
  rx <- matrix(c(1, 0, 0,
                 0, cos(el), -sin(el),
                 0, sin(el),  cos(el)), 3, 3, byrow = TRUE)
  rx %*% rz
}

make_parametric_mesh <- function(f, n_angle = 61L, n_vertical = 25L) {
  angle <- seq(0, 2 * pi, length.out = n_angle)
  vertical <- seq(-0.98, 0.98, length.out = n_vertical)
  grid <- expand.grid(angle = angle, vertical = vertical)
  template <- cbind(
    sqrt(1 - grid$vertical^2) * cos(grid$angle),
    sqrt(1 - grid$vertical^2) * sin(grid$angle),
    grid$vertical
  )
  list(values = f(template), n_angle = n_angle,
       n_vertical = n_vertical)
}

draw_projected_surface <- function(points, mesh, main, mesh_color,
                                   overlay = NULL,
                                   overlay_color = "#009E73") {
  rotation <- camera_rotation()
  all_values <- rbind(points, mesh$values, overlay)
  rotated_all <- t(rotation %*% t(all_values))
  limit <- max(abs(rotated_all[, 1:2])) * 1.05
  plot(NA, xlim = c(-limit, limit), ylim = c(-limit, limit),
       asp = 1, axes = FALSE, xlab = "", ylab = "", main = main)

  point_rotated <- t(rotation %*% t(points))
  order_points <- order(point_rotated[, 3])
  points(point_rotated[order_points, 1], point_rotated[order_points, 2],
         pch = 16, cex = 0.38, col = adjustcolor("grey35", 0.26))

  mesh_rotated <- t(rotation %*% t(mesh$values))
  mesh_matrix <- matrix(seq_len(nrow(mesh_rotated)),
                        nrow = mesh$n_angle, ncol = mesh$n_vertical)
  for (j in seq_len(mesh$n_vertical)) {
    ids <- mesh_matrix[, j]
    lines(mesh_rotated[ids, 1], mesh_rotated[ids, 2],
          col = adjustcolor(mesh_color, 0.75), lwd = 1.15)
  }
  meridians <- unique(round(seq(1, mesh$n_angle, length.out = 13L)))
  for (j in meridians) {
    ids <- mesh_matrix[j, ]
    lines(mesh_rotated[ids, 1], mesh_rotated[ids, 2],
          col = adjustcolor(mesh_color, 0.58), lwd = 0.85)
  }

  if (!is.null(overlay)) {
    overlay_rotated <- t(rotation %*% t(overlay))
    order_overlay <- order(overlay_rotated[, 3])
    points(overlay_rotated[order_overlay, 1],
           overlay_rotated[order_overlay, 2],
           pch = 16, cex = 0.55,
           col = adjustcolor(overlay_color, 0.78))
  }
  box()
}

true_mesh <- make_parametric_mesh(function(t) {
  flower2d3D_func(t, r0 = 1, a = 0.30, petals = 5, b = 0.50)
})
selected_mesh <- make_parametric_mesh(objects$full_pme$fit)

png(file.path(result_dir, "04-fixed-bandwidth-selected-surface.png"),
    width = 2600, height = 1200, res = 210)
par(mfrow = c(1, 2), mar = c(1, 1, 4.0, 1))
draw_projected_surface(
  objects$response_observed, true_mesh,
  "Observed responses and true flower", "black"
)
draw_projected_surface(
  objects$response_observed, selected_mesh,
  sprintf("Local-linear-selected CPME\nh = %.2f, lambda = %.2e",
          objects$bandwidth, objects$lambda),
  "#0072B2"
)
dev.off()

png(file.path(result_dir, "05-fixed-bandwidth-oof-predictions.png"),
    width = 2600, height = 1200, res = 210)
par(mfrow = c(1, 2), mar = c(1, 1, 4.0, 1))
draw_projected_surface(
  objects$response_observed, true_mesh,
  "Observed responses", "black"
)
draw_projected_surface(
  objects$response_observed, true_mesh,
  sprintf("Local-linear OOF predictions\nmean fold raw RMSE = %.4f",
          objects$mean_fold_raw_rmse),
  "black", overlay = objects$oof_prediction,
  overlay_color = "#009E73"
)
dev.off()

fixed <- sensitivity[abs(sensitivity$bandwidth - objects$bandwidth) < 1e-12, ]
methods <- c("kernel", "local_linear")
labels <- c(kernel = "Extrinsic kernel",
            local_linear = "Extrinsic local-linear")

png(file.path(result_dir, "06-fixed-bandwidth-lambda-performance.png"),
    width = 2800, height = 1250, res = 210)
par(mfrow = c(1, 2), mar = c(4.6, 4.7, 3.5, 1.0),
    mgp = c(2.7, 0.8, 0))
for (method in methods) {
  values <- fixed[fixed$method == method, ]
  values <- values[order(values$lambda), ]
  best <- values[which.min(values$mean_raw_rmse), ]
  y_range <- range(
    values$mean_raw_rmse + c(-1, 1) * values$se_raw_rmse,
    values$mean_oracle_rmse
  )
  plot(values$lambda, values$mean_raw_rmse,
       type = "b", log = "x", pch = 16, lwd = 2.2,
       col = "#D55E00", ylim = y_range,
       xlab = expression(lambda), ylab = "Five-fold RMSE",
       main = sprintf("%s, h = %.2f", labels[method], objects$bandwidth))
  arrows(values$lambda,
         values$mean_raw_rmse - values$se_raw_rmse,
         values$lambda,
         values$mean_raw_rmse + values$se_raw_rmse,
         angle = 90, code = 3, length = 0.035, col = "#D55E00")
  lines(values$lambda, values$mean_oracle_rmse,
        type = "b", pch = 1, lwd = 2.2, col = "#0072B2")
  abline(v = best$lambda, lty = 3, lwd = 1.5)
  legend("topleft",
         c("raw excluded responses", "noise-free oracle",
           sprintf("selected lambda = %.2e", best$lambda)),
         col = c("#D55E00", "#0072B2", "black"),
         pch = c(16, 1, NA), lty = c(1, 1, 3),
         lwd = c(2.2, 2.2, 1.5), bty = "n", cex = 0.80)
}
dev.off()

cat(sprintf(
  paste0("Rendered h=%.2f, lambda=%.8e, mean fold raw RMSE=%.5f, ",
         "pooled raw RMSE=%.5f\n"),
  objects$bandwidth, objects$lambda, objects$mean_fold_raw_rmse,
  objects$raw_rmse
))
