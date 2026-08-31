## Hold the extrinsic-regression bandwidth fixed and determine which CPME
## lambda minimizes five-fold raw-response prediction error. Uses cached PME
## fits from generate-results.R and does not refit any manifolds.

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
cache_dir <- file.path(notebook_dir, "cache")

local_library <- file.path(repo_root, "packages")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
source(file.path(repo_root, "CompactPME", "R", "qm.R"))
source(file.path(repo_root, "CompactPME", "R", "project_optimize2.R"))

objects <- readRDS(file.path(
  result_dir, "extrinsic-regression-2d-flower-objects.rds"
))
z_external <- objects$data$z_external
response_observed <- objects$data$response_observed
conditional_mean <- objects$data$conditional_mean
fold_id <- objects$data$fold_id
lambda_grid <- objects$tuning$lambda_grid
bandwidth_grid <- objects$tuning$bandwidth_grid
methods <- c("kernel", "local_linear")
surface_grid <- fibonaccisphere(1800L)

gaussian_weights <- function(z_train, z_new, bandwidth) {
  squared_distance <- outer(rowSums(z_new^2), rowSums(z_train^2), "+") -
    2 * z_new %*% t(z_train)
  exp(-0.5 * pmax(squared_distance, 0) / bandwidth^2)
}

ambient_kernel <- function(z_train, y_train, z_new, bandwidth) {
  weights <- gaussian_weights(z_train, z_new, bandwidth)
  weights %*% y_train / rowSums(weights)
}

ambient_local_linear <- function(z_train, y_train, z_new, bandwidth) {
  output <- matrix(NA_real_, nrow = nrow(z_new), ncol = ncol(y_train))
  for (j in seq_len(nrow(z_new))) {
    delta <- sweep(z_train, 2L, z_new[j, ], "-")
    weights <- exp(-0.5 * rowSums(delta^2) / bandwidth^2)
    design <- cbind(1, delta)
    normal_matrix <- crossprod(design, weights * design)
    right_hand_side <- crossprod(design, weights * y_train)
    scale_value <- max(1, max(abs(diag(normal_matrix))))
    coefficient <- if (is.finite(rcond(normal_matrix)) &&
                       rcond(normal_matrix) > 1e-10) {
      tryCatch(solve(normal_matrix, right_hand_side),
               error = function(e) NULL)
    } else {
      NULL
    }
    if (is.null(coefficient)) {
      coefficient <- tryCatch(
        solve(normal_matrix + diag(1e-8 * scale_value, 3L),
              right_hand_side),
        error = function(e) NULL
      )
    }
    if (is.null(coefficient) || any(!is.finite(coefficient))) {
      output[j, ] <- colSums(weights * y_train) / sum(weights)
    } else {
      output[j, ] <- coefficient[1, ]
    }
  }
  output
}

project_grid <- function(points, surface) {
  surface_norm <- rowSums(surface^2)
  distance_squared <- outer(rowSums(points^2), surface_norm, "+") -
    2 * points %*% t(surface)
  surface[max.col(-distance_squared, ties.method = "first"), , drop = FALSE]
}

rmse_3d <- function(observed, predicted) {
  sqrt(mean(rowSums((observed - predicted)^2)))
}

rows <- vector(
  "list",
  length(lambda_grid) * length(bandwidth_grid) * length(methods) * 5L
)
row_id <- 0L

for (li in seq_along(lambda_grid)) {
  for (fold in 1:5) {
    training_id <- which(fold_id != fold)
    validation_id <- which(fold_id == fold)
    pme <- readRDS(file.path(
      cache_dir, sprintf("pme-lambda-%02d-fold-%d.rds", li, fold)
    ))
    surface <- pme$fit(surface_grid)

    for (bandwidth in bandwidth_grid) {
      for (method in methods) {
        ambient <- switch(
          method,
          kernel = ambient_kernel(
            z_external[training_id, , drop = FALSE], pme$denoised,
            z_external[validation_id, , drop = FALSE], bandwidth
          ),
          local_linear = ambient_local_linear(
            z_external[training_id, , drop = FALSE], pme$denoised,
            z_external[validation_id, , drop = FALSE], bandwidth
          )
        )
        prediction <- project_grid(ambient, surface)
        row_id <- row_id + 1L
        rows[[row_id]] <- data.frame(
          method = method,
          lambda = lambda_grid[li],
          bandwidth = bandwidth,
          fold = fold,
          raw_rmse = rmse_3d(
            response_observed[validation_id, , drop = FALSE], prediction
          ),
          oracle_rmse = rmse_3d(
            conditional_mean[validation_id, , drop = FALSE], prediction
          )
        )
      }
    }
  }
  cat(sprintf("Completed lambda %.3e\n", lambda_grid[li]))
}

fold_sensitivity <- do.call(rbind, rows)
groups <- split(
  fold_sensitivity,
  interaction(
    fold_sensitivity$method,
    fold_sensitivity$lambda,
    fold_sensitivity$bandwidth,
    drop = TRUE
  )
)
aggregate_sensitivity <- do.call(rbind, lapply(groups, function(group) {
  data.frame(
    method = group$method[1],
    lambda = group$lambda[1],
    bandwidth = group$bandwidth[1],
    mean_raw_rmse = mean(group$raw_rmse),
    se_raw_rmse = sd(group$raw_rmse) / sqrt(nrow(group)),
    mean_oracle_rmse = mean(group$oracle_rmse)
  )
}))
rownames(aggregate_sensitivity) <- NULL
aggregate_sensitivity <- aggregate_sensitivity[
  order(
    aggregate_sensitivity$method,
    aggregate_sensitivity$bandwidth,
    aggregate_sensitivity$lambda
  ),
]

bandwidth_groups <- split(
  aggregate_sensitivity,
  interaction(
    aggregate_sensitivity$method,
    aggregate_sensitivity$bandwidth,
    drop = TRUE
  )
)
selected_by_bandwidth <- do.call(rbind, lapply(bandwidth_groups, function(x) {
  x[which.min(x$mean_raw_rmse), , drop = FALSE]
}))
rownames(selected_by_bandwidth) <- NULL
selected_by_bandwidth <- selected_by_bandwidth[
  order(selected_by_bandwidth$method, selected_by_bandwidth$bandwidth),
]

write.csv(
  fold_sensitivity,
  file.path(result_dir, "bandwidth-sensitivity-fold-performance.csv"),
  row.names = FALSE
)
write.csv(
  aggregate_sensitivity,
  file.path(result_dir, "bandwidth-sensitivity-all-combinations.csv"),
  row.names = FALSE
)
write.csv(
  selected_by_bandwidth,
  file.path(result_dir, "lambda-selected-by-fixed-bandwidth.csv"),
  row.names = FALSE
)

print(selected_by_bandwidth)
