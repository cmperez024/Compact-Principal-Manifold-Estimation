# Migrate legacy cached RMSE summaries to the draft's mean fold-MSE criterion.
# Squaring the mean RMSE is not the same as averaging the squared fold RMSEs.
upgrade_cv_mse <- function(object) {
  if (is.null(object$tuning)) return(object)
  tuning <- object$tuning
  fold_name <- if (!is.null(tuning$fold_scores)) "fold_scores" else "fold_performance"
  summary_name <- if (!is.null(tuning$aggregate_scores)) "aggregate_scores" else "lambda_performance"
  folds <- tuning[[fold_name]]
  scores <- tuning[[summary_name]]
  if (is.null(folds) || is.null(scores)) return(object)
  for (i in seq_len(nrow(scores))) {
    ids <- folds$method == scores$method[i] & folds$lambda == scores$lambda[i]
    raw <- folds$validation_rmse_raw[ids]^2
    oracle <- folds$validation_rmse_oracle[ids]^2
    if (!length(raw)) stop("Missing fold errors for a cached CV summary.")
    scores$mean_raw_mse[i] <- mean(raw)
    scores$se_raw_mse[i] <- stats::sd(raw)/sqrt(length(raw))
    scores$mean_oracle_mse[i] <- mean(oracle)
    scores$se_oracle_mse[i] <- stats::sd(oracle)/sqrt(length(oracle))
  }
  scores$selected <- FALSE
  for (method in unique(scores$method)) {
    ids <- which(scores$method == method)
    best <- ids[which.min(scores$mean_raw_mse[ids])]
    choice <- scores[best, , drop = FALSE]
    previous <- object$final[[method]]$choice$lambda
    if (!is.null(previous) && !isTRUE(all.equal(previous, choice$lambda, tolerance = 0))) {
      stop("MSE selects a different lambda from this legacy cached final fit. ",
           "Regenerate the experiment in a new PME_RESULTS_ROOT directory.")
    }
    scores$selected[best] <- TRUE
    tuning$selected[[method]] <- scores[best, , drop = FALSE]
    if (!is.null(object$final[[method]])) {
      object$final[[method]]$choice <- scores[best, , drop = FALSE]
    }
  }
  tuning[[summary_name]] <- scores
  tuning$selection_loss <- "mean_fold_mse"
  object$tuning <- tuning
  object
}

# Resolve paths from this file, including when invoked outside the project.
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_file)) stop("Run this script with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_file[1L]),
                                     winslash = "/", mustWork = TRUE))

## Static figures for h = 0.06 and the local-linear selected lambda.

options(stringsAsFactors = FALSE)

find_repo_root <- function(start = script_dir) {
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
notebook_dir <- Sys.getenv("PME_RESULTS_ROOT", unset = notebook_dir)
result_dir <- file.path(notebook_dir, "results")
local_library <- file.path(repo_root, "packages")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))

cart2sph <- pracma::cart2sph
sph2cart <- pracma::sph2cart
if (!requireNamespace("CompactPME", quietly = TRUE)) {
  stop("Install the CompactPME package before running this script.")
}
spline2d <- CompactPME::spline2d
qm <- CompactPME:::qm
expr <- CompactPME:::expr
krt <- CompactPME:::krt
krst <- CompactPME:::krst
Rk <- CompactPME:::Rk
project_optimize2 <- CompactPME:::project_optimize2
fibonaccisphere <- CompactPME:::fibonaccisphere
flower2d3D_func <- function(t, r0 = 1, a = 0.3, petals = 6, b = 0.5) {
  t <- as.matrix(t)
  t <- t / sqrt(rowSums(t^2))
  longitude <- atan2(t[, 2L], t[, 1L])
  z <- pmax(-1, pmin(1, t[, 3L]))
  radius_xy <- r0 * (1 + a * cos(petals * longitude))
  cbind(radius_xy * sqrt(1 - z^2) * cos(longitude),
        radius_xy * sqrt(1 - z^2) * sin(longitude), b * z)
}
# This script includes its own data/PME generation fallback.
reproduce_inputs <- function(target_dir, variant = "original", n_sample = 400L) {
options(stringsAsFactors = FALSE)
set.seed(20260807)

arguments <- c(variant, as.character(n_sample))
experiment_variant <- if (length(arguments)) arguments[[1L]] else "original"
if (!experiment_variant %in% c("original", "uniform-noninjective")) {
  stop("Variant must be 'original' or 'uniform-noninjective'.")
}
sample_size <- if (length(arguments) >= 2L) {
  as.integer(arguments[[2L]])
} else {
  400L
}
if (!is.finite(sample_size) || sample_size < 10L || sample_size %% 5L != 0L) {
  stop("Sample size must be an integer of at least 10 and divisible by 5.")
}

notebook_dir <- target_dir
result_suffix <- if (experiment_variant == "original") {
  ""
} else {
  sprintf("-uniform-n%d", sample_size)
}
cache_suffix <- if (experiment_variant == "original") {
  ""
} else {
  sprintf("-uniform-noninjective-n%d", sample_size)
}
notebook_dir <- Sys.getenv("PME_RESULTS_ROOT", unset = notebook_dir)
if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) sample_size <- 30L
result_dir <- file.path(notebook_dir, paste0("results", result_suffix))
cache_dir <- file.path(notebook_dir, paste0("cache", cache_suffix))
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

local_library <- file.path(repo_root, "packages")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
if (!requireNamespace("CompactPME", quietly = TRUE)) {
  stop("Install the CompactPME package before running this script.")
}

## All project functionality comes from the installed CompactPME package.
spline2d <- CompactPME::spline2d
project_optimize2 <- CompactPME:::project_optimize2
fibonaccisphere <- CompactPME:::fibonaccisphere

flower2d3D_func <- function(t, r0 = 1, a = 0.3, petals = 6, b = 0.5) {
  t <- as.matrix(t)
  t <- t / sqrt(rowSums(t^2))
  longitude <- atan2(t[, 2L], t[, 1L])
  z <- pmax(-1, pmin(1, t[, 3L]))
  radius_xy <- r0 * (1 + a * cos(petals * longitude))
  cbind(
    radius_xy * sqrt(1 - z^2) * cos(longitude),
    radius_xy * sqrt(1 - z^2) * sin(longitude),
    b * z
  )
}

normalize_rows <- function(x) {
  x <- as.matrix(x)
  norms <- sqrt(rowSums(x^2))
  if (any(norms <= .Machine$double.eps)) {
    stop("Cannot normalize a zero-length row.")
  }
  x / norms
}

phase_map_s2 <- function(z) {
  z <- as.matrix(z)
  if (experiment_variant == "uniform-noninjective") {
    u1 <- (z[, 1] + 1) / 2
    u2 <- (z[, 2] + 1) / 2
    q <- (2 * u1 + 0.22 * sin(2 * pi * u2)) %% 1
    vertical <- (2 / pi) * asin(sin(2 * pi * u2))
  } else {
    q <- (0.5 + 0.45 * sin(
      pi * z[, 1] + 0.35 * sin(pi * z[, 2])
    )) %% 1
    vertical <- 0.85 * sin(
      pi * z[, 2] + 0.25 * sin(pi * z[, 1])
    )
  }
  angle <- 2 * pi * q
  cbind(
    sqrt(1 - vertical^2) * cos(angle),
    sqrt(1 - vertical^2) * sin(angle),
    vertical
  )
}

fit_pme_s2 <- function(x, lambda, min_iterations = 3L,
                       max_iterations = if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 2L else 15L,
                       relative_tolerance = 1e-4) {
  centered <- scale(x, scale = FALSE)
  projection <- normalize_rows(centered)
  fit <- spline2d(x, projection, lambda)
  projection <- project_optimize2(x, fit, projection)
  fidelity <- mean(rowSums((x - fit(projection))^2))
  history <- data.frame(
    iteration = 1L,
    fidelity = fidelity,
    relative_change = NA_real_
  )
  converged <- FALSE

  for (iteration in 2:max_iterations) {
    next_fit <- spline2d(x, projection, lambda)
    next_projection <- project_optimize2(x, next_fit, projection)
    next_fidelity <- mean(rowSums((x - next_fit(next_projection))^2))
    relative_change <- abs(next_fidelity - fidelity) / max(fidelity, 1e-14)
    history <- rbind(history, data.frame(
      iteration = iteration,
      fidelity = next_fidelity,
      relative_change = relative_change
    ))
    fit <- next_fit
    projection <- next_projection
    fidelity <- next_fidelity
    if (iteration >= min_iterations &&
        relative_change <= relative_tolerance) {
      converged <- TRUE
      break
    }
  }

  list(
    fit = fit,
    projection_index = projection,
    denoised = fit(projection),
    history = history,
    converged = converged,
    lambda = lambda
  )
}

surface_grid <- fibonaccisphere(if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 100L else 1800L)

make_surface_cache <- function(fit) {
  list(template = surface_grid, ambient = fit(surface_grid))
}

project_to_surface_grid <- function(points, surface_cache,
                                    block_size = 250L) {
  points <- as.matrix(points)
  surface <- surface_cache$ambient
  nearest <- integer(nrow(points))
  surface_norm <- rowSums(surface^2)
  starts <- seq.int(1L, nrow(points), by = block_size)
  for (start in starts) {
    stop_at <- min(nrow(points), start + block_size - 1L)
    ids <- start:stop_at
    distance_squared <- outer(rowSums(points[ids, , drop = FALSE]^2),
                              surface_norm, "+") -
      2 * points[ids, , drop = FALSE] %*% t(surface)
    nearest[ids] <- max.col(-distance_squared, ties.method = "first")
  }
  list(
    index = surface_cache$template[nearest, , drop = FALSE],
    fitted = surface[nearest, , drop = FALSE]
  )
}

gaussian_weights <- function(z_train, z_new, bandwidth) {
  z_train <- as.matrix(z_train)
  z_new <- as.matrix(z_new)
  squared_distance <- outer(rowSums(z_new^2), rowSums(z_train^2), "+") -
    2 * z_new %*% t(z_train)
  exp(-0.5 * pmax(squared_distance, 0) / bandwidth^2)
}

ambient_kernel <- function(z_train, y_train, z_new, bandwidth,
                           leave_one_out = FALSE) {
  weights <- gaussian_weights(z_train, z_new, bandwidth)
  if (leave_one_out) {
    if (nrow(z_train) != nrow(z_new) ||
        !isTRUE(all.equal(z_train, z_new))) {
      stop("LOO kernel prediction requires identical training and new Z.")
    }
    diag(weights) <- 0
  }
  weight_sums <- rowSums(weights)
  if (any(weight_sums <= .Machine$double.eps)) {
    stop("Bandwidth produced effectively zero total kernel weight.")
  }
  weights %*% y_train / weight_sums
}

ambient_local_linear <- function(z_train, y_train, z_new, bandwidth,
                                 leave_one_out = FALSE) {
  z_train <- as.matrix(z_train)
  z_new <- as.matrix(z_new)
  y_train <- as.matrix(y_train)
  output <- matrix(NA_real_, nrow = nrow(z_new), ncol = ncol(y_train))

  for (j in seq_len(nrow(z_new))) {
    delta <- sweep(z_train, 2L, z_new[j, ], "-")
    weights <- exp(-0.5 * rowSums(delta^2) / bandwidth^2)
    if (leave_one_out) weights[j] <- 0
    design <- cbind(1, delta)
    normal_matrix <- crossprod(design, weights * design)
    right_hand_side <- crossprod(design, weights * y_train)
    scale_value <- max(1, max(abs(diag(normal_matrix))))
    stable <- is.finite(rcond(normal_matrix)) && rcond(normal_matrix) > 1e-10
    coefficient <- if (stable) {
      tryCatch(
        solve(normal_matrix, right_hand_side),
        error = function(e) NULL
      )
    } else {
      NULL
    }
    if (is.null(coefficient)) {
      regularized <- normal_matrix + diag(1e-8 * scale_value, 3L)
      coefficient <- tryCatch(
        solve(regularized, right_hand_side),
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

extrinsic_predict <- function(method, z_train, y_train, z_new, bandwidth,
                              surface_cache, leave_one_out = FALSE) {
  ambient <- switch(
    method,
    kernel = ambient_kernel(
      z_train, y_train, z_new, bandwidth, leave_one_out
    ),
    local_linear = ambient_local_linear(
      z_train, y_train, z_new, bandwidth, leave_one_out
    ),
    stop("Unknown regression method.")
  )
  projection <- project_to_surface_grid(ambient, surface_cache)
  list(ambient = ambient, index = projection$index,
       extrinsic = projection$fitted)
}

rmse_3d <- function(observed, predicted) {
  sqrt(mean(rowSums((observed - predicted)^2)))
}

select_bandwidth <- function(method, z_train, denoised_train,
                             bandwidth_grid, surface_cache) {
  error <- vapply(bandwidth_grid, function(bandwidth) {
    prediction <- extrinsic_predict(
      method, z_train, denoised_train, z_train, bandwidth,
      surface_cache, leave_one_out = TRUE
    )$extrinsic
    rmse_3d(denoised_train, prediction)
  }, numeric(1))
  selected <- which.min(error)
  list(
    bandwidth = bandwidth_grid[selected],
    rmse = error[selected],
    errors = error
  )
}

## Fixed-design response-noise experiment -----------------------------------
n <- sample_size
k_folds <- 5L
noise_sd <- 0.04
z_external <- cbind(
  z1 = runif(n, -1, 1),
  z2 = runif(n, -1, 1)
)
latent_template <- phase_map_s2(z_external)
conditional_mean <- flower2d3D_func(
  latent_template, r0 = 1, a = 0.30, petals = 5, b = 0.50
)
response_error <- matrix(rnorm(3L * n, sd = noise_sd), nrow = n, ncol = 3L)
response_observed <- conditional_mean + response_error
fold_id <- sample(rep(seq_len(k_folds), length.out = n))
stopifnot(all(tabulate(fold_id, nbins = k_folds) == n / k_folds))

lambda_grid <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) c(1e-6, 1e-4) else c(10^seq(-12, -4, length.out = 8L), 1e-3, 1e-2)
bandwidth_grid <- c(0.03, 0.04, 0.05, 0.06, 0.08, 0.12,
                    0.18, 0.26, 0.38, 0.55, 0.80)
if (identical(Sys.getenv("PME_SMOKE_TEST"), "1") && exists("bandwidth_grid")) {
  bandwidth_grid <- c(0.2, 0.4, 0.8)
}
methods <- c("kernel", "local_linear")
method_labels <- c(
  kernel = "Extrinsic kernel",
  local_linear = "Extrinsic local-linear"
)

prediction_store <- array(
  NA_real_,
  dim = c(n, 3L, length(lambda_grid), length(methods)),
  dimnames = list(NULL, c("X1", "X2", "X3"),
                  format(lambda_grid, scientific = TRUE), methods)
)
fold_rows <- vector("list", length(lambda_grid) * k_folds * length(methods))
fit_rows <- vector("list", length(lambda_grid) * k_folds)
row_id <- 0L
fit_id <- 0L

for (li in seq_along(lambda_grid)) {
  lambda <- lambda_grid[li]
  for (fold in seq_len(k_folds)) {
    validation_id <- which(fold_id == fold)
    training_id <- which(fold_id != fold)
    cache_path <- file.path(
      cache_dir, sprintf("pme-lambda-%02d-fold-%d.rds", li, fold)
    )
    pme <- if (file.exists(cache_path)) {
      readRDS(cache_path)
    } else {
      value <- fit_pme_s2(
        response_observed[training_id, , drop = FALSE], lambda
      )
      saveRDS(value, cache_path)
      value
    }
    fit_id <- fit_id + 1L
    fit_rows[[fit_id]] <- data.frame(
      lambda = lambda,
      fold = fold,
      iterations = nrow(pme$history),
      converged = pme$converged,
      final_fidelity = tail(pme$history$fidelity, 1L)
    )
    surface_cache_fold <- make_surface_cache(pme$fit)

    for (mi in seq_along(methods)) {
      method <- methods[mi]
      bandwidth_selection <- select_bandwidth(
        method,
        z_external[training_id, , drop = FALSE],
        pme$denoised,
        bandwidth_grid,
        surface_cache_fold
      )
      prediction <- extrinsic_predict(
        method,
        z_external[training_id, , drop = FALSE],
        pme$denoised,
        z_external[validation_id, , drop = FALSE],
        bandwidth_selection$bandwidth,
        surface_cache_fold
      )$extrinsic
      prediction_store[validation_id, , li, mi] <- prediction
      row_id <- row_id + 1L
      fold_rows[[row_id]] <- data.frame(
        method = method,
        lambda = lambda,
        fold = fold,
        bandwidth = bandwidth_selection$bandwidth,
        training_loo_rmse = bandwidth_selection$rmse,
        validation_rmse_raw = rmse_3d(
          response_observed[validation_id, , drop = FALSE], prediction
        ),
        validation_rmse_oracle = rmse_3d(
          conditional_mean[validation_id, , drop = FALSE], prediction
        )
      )
    }
    cat(sprintf(
      "lambda=%9.2e fold=%d iterations=%d converged=%s\n",
      lambda, fold, nrow(pme$history), pme$converged
    ))
  }
}

fold_performance <- do.call(rbind, fold_rows)
pme_status <- do.call(rbind, fit_rows)

aggregate_one <- function(data) {
  raw <- data$validation_rmse_raw
  oracle <- data$validation_rmse_oracle
  data.frame(
    mean_raw_mse = mean(raw^2),
    se_raw_mse = sd(raw^2) / sqrt(length(raw)),
    mean_oracle_mse = mean(oracle^2),
    se_oracle_mse = sd(oracle^2) / sqrt(length(oracle)),
    mean_raw_rmse = mean(raw),
    se_raw_rmse = sd(raw) / sqrt(length(raw)),
    mean_oracle_rmse = mean(oracle),
    se_oracle_rmse = sd(oracle) / sqrt(length(oracle)),
    mean_bandwidth = mean(data$bandwidth)
  )
}

groups <- split(
  fold_performance,
  interaction(fold_performance$method, fold_performance$lambda, drop = TRUE)
)
lambda_performance <- do.call(rbind, lapply(groups, function(group) {
  cbind(
    data.frame(method = group$method[1], lambda = group$lambda[1]),
    aggregate_one(group)
  )
}))
rownames(lambda_performance) <- NULL
lambda_performance <- lambda_performance[
  order(lambda_performance$method, lambda_performance$lambda),
]
lambda_performance$selected <- FALSE

selected <- setNames(vector("list", length(methods)), methods)
for (method in methods) {
  candidates <- which(lambda_performance$method == method)
  best <- candidates[which.min(lambda_performance$mean_raw_mse[candidates])]
  lambda_performance$selected[best] <- TRUE
  selected[[method]] <- lambda_performance[best, , drop = FALSE]
}

## Full-data refits and bandwidth choices -----------------------------------
final <- setNames(vector("list", length(methods)), methods)
full_fit_by_lambda <- list()
for (method in methods) {
  choice <- selected[[method]]
  lambda <- choice$lambda
  cache_key <- format(lambda, scientific = TRUE, digits = 16)
  if (is.null(full_fit_by_lambda[[cache_key]])) {
    full_fit_by_lambda[[cache_key]] <- fit_pme_s2(response_observed, lambda)
  }
  pme <- full_fit_by_lambda[[cache_key]]
  surface_cache_full <- make_surface_cache(pme$fit)
  bandwidth_selection <- select_bandwidth(
    method, z_external, pme$denoised, bandwidth_grid, surface_cache_full
  )
  final[[method]] <- list(
    choice = choice,
    bandwidth = bandwidth_selection$bandwidth,
    bandwidth_selection = bandwidth_selection,
    pme = pme,
    surface_cache = surface_cache_full,
    oof_prediction = prediction_store[, ,
      which(lambda_grid == lambda), which(methods == method)]
  )
}

write.csv(
  fold_performance,
  file.path(result_dir, "fold-performance.csv"),
  row.names = FALSE
)
write.csv(
  lambda_performance,
  file.path(result_dir, "lambda-performance.csv"),
  row.names = FALSE
)
write.csv(
  pme_status,
  file.path(result_dir, "pme-fit-status.csv"),
  row.names = FALSE
)
saveRDS(list(
  model = list(
    n = n,
    k_folds = k_folds,
    noise_sd = noise_sd,
    flower = list(r0 = 1, a = 0.30, petals = 5L, b = 0.50),
    phase_map = if (experiment_variant == "uniform-noninjective") {
      paste(
        "u=(z+1)/2; q=(2*u1+0.22*sin(2*pi*u2)) mod 1;",
        "s=(2/pi)*asin(sin(2*pi*u2))"
      )
    } else {
      paste(
        "q=(0.5+0.45*sin(pi*z1+0.35*sin(pi*z2))) mod 1;",
        "s=0.85*sin(pi*z2+0.25*sin(pi*z1))"
      )
    },
    experiment_variant = experiment_variant
  ),
  data = list(
    z_external = z_external,
    latent_template = latent_template,
    conditional_mean = conditional_mean,
    response_error = response_error,
    response_observed = response_observed,
    fold_id = fold_id
  ),
  tuning = list(
    lambda_grid = lambda_grid,
    bandwidth_grid = bandwidth_grid,
    fold_performance = fold_performance,
    lambda_performance = lambda_performance,
    selected = selected
  ),
  predictions = prediction_store,
  final = final
), file.path(result_dir, "extrinsic-regression-2d-flower-objects.rds"))

## Figure 1: downstream CV performance --------------------------------------
png(
  file.path(result_dir, "01-kfold-lambda-performance.png"),
  width = 2800, height = 1250, res = 210
)
par(mfrow = c(1, 2), mar = c(4.5, 4.7, 3.4, 1.0),
    mgp = c(2.7, 0.8, 0))
for (method in methods) {
  values <- lambda_performance[lambda_performance$method == method, ]
  y_range <- range(
    values$mean_raw_mse + c(-1, 1) * values$se_raw_mse,
    values$mean_oracle_mse + c(-1, 1) * values$se_oracle_mse
  )
  plot(
    values$lambda, values$mean_raw_mse,
    log = "x", type = "b", pch = 16, lwd = 2.2,
    col = "#D55E00", ylim = y_range,
    xlab = expression(lambda), ylab = "Five-fold MSE",
    main = method_labels[method]
  )
  arrows(
    values$lambda, values$mean_raw_mse - values$se_raw_mse,
    values$lambda, values$mean_raw_mse + values$se_raw_mse,
    angle = 90, code = 3, length = 0.035, col = "#D55E00"
  )
  lines(values$lambda, values$mean_oracle_mse,
        type = "b", pch = 1, lwd = 2.2, col = "#0072B2")
  arrows(
    values$lambda, values$mean_oracle_mse - values$se_oracle_mse,
    values$lambda, values$mean_oracle_mse + values$se_oracle_mse,
    angle = 90, code = 3, length = 0.035, col = "#0072B2"
  )
  best <- values[values$selected, ]
  abline(v = best$lambda, lty = 3, lwd = 1.5)
  legend(
    "topleft",
    c("raw excluded responses", "noise-free oracle", "selected lambda"),
    col = c("#D55E00", "#0072B2", "black"),
    pch = c(16, 1, NA), lty = c(1, 1, 3), lwd = c(2.2, 2.2, 1.5),
    bty = "n", cex = 0.82
  )
}
dev.off()

## Helpers for static 3D comparison figures --------------------------------
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
  list(
    values = f(template),
    n_angle = n_angle,
    n_vertical = n_vertical
  )
}

draw_projected_surface <- function(points, mesh, main, mesh_color,
                                   overlay = NULL, overlay_color = NULL) {
  rotation <- camera_rotation()
  all_values <- rbind(points, mesh$values, overlay)
  rotated_all <- t(rotation %*% t(all_values))
  limits <- range(rotated_all[, 1:2]) * 1.04
  plot(NA, xlim = limits, ylim = limits, asp = 1, axes = FALSE,
       xlab = "", ylab = "", main = main)

  mesh_rotated <- t(rotation %*% t(mesh$values))
  point_rotated <- t(rotation %*% t(points))
  point_order <- order(point_rotated[, 3])
  points(point_rotated[point_order, 1], point_rotated[point_order, 2],
         pch = 16, cex = 0.36, col = adjustcolor("grey35", 0.28))

  mesh_matrix <- matrix(
    seq_len(nrow(mesh_rotated)),
    nrow = mesh$n_angle,
    ncol = mesh$n_vertical
  )
  for (j in seq_len(mesh$n_vertical)) {
    ids <- mesh_matrix[, j]
    lines(mesh_rotated[ids, 1], mesh_rotated[ids, 2],
          col = adjustcolor(mesh_color, 0.72), lwd = 1.1)
  }
  meridians <- unique(round(seq(1, mesh$n_angle, length.out = 13L)))
  for (j in meridians) {
    ids <- mesh_matrix[j, ]
    lines(mesh_rotated[ids, 1], mesh_rotated[ids, 2],
          col = adjustcolor(mesh_color, 0.55), lwd = 0.8)
  }
  if (!is.null(overlay)) {
    overlay_rotated <- t(rotation %*% t(overlay))
    overlay_order <- order(overlay_rotated[, 3])
    points(overlay_rotated[overlay_order, 1], overlay_rotated[overlay_order, 2],
           pch = 16, cex = 0.45,
           col = adjustcolor(overlay_color, 0.72))
  }
  box()
}

true_mesh <- make_parametric_mesh(function(t) {
  flower2d3D_func(t, r0 = 1, a = 0.30, petals = 5, b = 0.50)
})

png(
  file.path(result_dir, "02-data-and-selected-pme-surfaces.png"),
  width = 3200, height = 1150, res = 210
)
par(mfrow = c(1, 3), mar = c(1.0, 1.0, 3.7, 1.0))
draw_projected_surface(
  response_observed, true_mesh,
  "Observed data and true surface", "black"
)
surface_colors <- c(kernel = "#0072B2", local_linear = "#009E73")
for (method in methods) {
  selected_mesh <- make_parametric_mesh(final[[method]]$pme$fit)
  draw_projected_surface(
    response_observed, selected_mesh,
    sprintf(
      "%s-selected PME\nlambda = %.2e",
      method_labels[method], final[[method]]$choice$lambda
    ),
    surface_colors[method]
  )
}
dev.off()

png(
  file.path(result_dir, "03-out-of-fold-predictions.png"),
  width = 2500, height = 1200, res = 210
)
par(mfrow = c(1, 2), mar = c(1.0, 1.0, 3.7, 1.0))
for (method in methods) {
  draw_projected_surface(
    response_observed, true_mesh,
    sprintf(
      "%s OOF predictions\nraw RMSE = %.4f",
      method_labels[method], final[[method]]$choice$mean_raw_rmse
    ),
    "black",
    overlay = final[[method]]$oof_prediction,
    overlay_color = surface_colors[method]
  )
}
dev.off()

for (method in methods) {
  choice <- selected[[method]]
  cat(sprintf(
    "%s: lambda=%.3e, mean OOF raw RMSE=%.5f, oracle RMSE=%.5f, full h=%.2f\n",
    method_labels[method], choice$lambda, choice$mean_raw_rmse,
    choice$mean_oracle_rmse, final[[method]]$bandwidth
  ))
}

}
fibonaccisphere <- CompactPME:::fibonaccisphere
qm <- CompactPME:::qm
expr <- CompactPME:::expr
input_file <- file.path(result_dir, "extrinsic-regression-2d-flower-objects.rds")
if (!file.exists(input_file)) reproduce_inputs(notebook_dir, variant = "original", n_sample = 400L)

cache_dir <- file.path(notebook_dir, "cache")
reproduce_fixed_bandwidth <- function() {
objects <- readRDS(file.path(
  result_dir, "extrinsic-regression-2d-flower-objects.rds"
))
objects <- upgrade_cv_mse(objects)
z_external <- objects$data$z_external
response_observed <- objects$data$response_observed
fold_id <- objects$data$fold_id
lambda_grid <- objects$tuning$lambda_grid
target_lambda <- lambda_grid[if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) length(lambda_grid) else 7L]
target_bandwidth <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 0.2 else 0.06
lambda_index <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) length(lambda_grid) else 7L

normalize_rows <- function(x) x / sqrt(rowSums(x^2))

fit_pme_s2 <- function(x, lambda, max_iterations = if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 2L else 15L) {
  projection <- normalize_rows(scale(x, scale = FALSE))
  fit <- spline2d(x, projection, lambda)
  projection <- project_optimize2(x, fit, projection)
  for (iteration in 2:max_iterations) {
    fit <- spline2d(x, projection, lambda)
    projection <- project_optimize2(x, fit, projection)
  }
  list(fit = fit, projection = projection, denoised = fit(projection))
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

surface_grid <- fibonaccisphere(if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 100L else 1800L)
project_grid <- function(points, surface) {
  distance_squared <- outer(rowSums(points^2), rowSums(surface^2), "+") -
    2 * points %*% t(surface)
  surface[max.col(-distance_squared, ties.method = "first"), , drop = FALSE]
}

## Honest OOF predictions with h fixed at 0.06.
oof_prediction <- matrix(NA_real_, nrow(response_observed), 3L)
for (fold in 1:5) {
  training_id <- which(fold_id != fold)
  validation_id <- which(fold_id == fold)
  pme <- readRDS(file.path(
    cache_dir,
    sprintf("pme-lambda-%02d-fold-%d.rds", lambda_index, fold)
  ))
  ambient <- ambient_local_linear(
    z_external[training_id, , drop = FALSE],
    pme$denoised,
    z_external[validation_id, , drop = FALSE],
    target_bandwidth
  )
  oof_prediction[validation_id, ] <- project_grid(
    ambient, pme$fit(surface_grid)
  )
}

## Full-data surface at the same lambda for geometric display only.
full_pme <- fit_pme_s2(response_observed, target_lambda)

fold_raw_rmse <- vapply(1:5, function(fold) {
  ids <- which(fold_id == fold)
  sqrt(mean(rowSums(
    (response_observed[ids, , drop = FALSE] -
       oof_prediction[ids, , drop = FALSE])^2
  )))
}, numeric(1))
fold_oracle_rmse <- vapply(1:5, function(fold) {
  ids <- which(fold_id == fold)
  sqrt(mean(rowSums(
    (objects$data$conditional_mean[ids, , drop = FALSE] -
       oof_prediction[ids, , drop = FALSE])^2
  )))
}, numeric(1))

saveRDS(list(
  lambda = target_lambda,
  bandwidth = target_bandwidth,
  z_external = z_external,
  response_observed = response_observed,
  conditional_mean = objects$data$conditional_mean,
  fold_id = fold_id,
  oof_prediction = oof_prediction,
  raw_rmse = sqrt(mean(rowSums((response_observed - oof_prediction)^2))),
  oracle_rmse = sqrt(mean(rowSums(
    (objects$data$conditional_mean - oof_prediction)^2
  ))),
  fold_raw_rmse = fold_raw_rmse,
  fold_oracle_rmse = fold_oracle_rmse,
  mean_fold_raw_rmse = mean(fold_raw_rmse),
  mean_fold_oracle_rmse = mean(fold_oracle_rmse),
  full_pme = full_pme
), file.path(result_dir, "fixed-bandwidth-local-linear-objects.rds"))

}
reproduce_sensitivity <- function() {
objects <- readRDS(file.path(
  result_dir, "extrinsic-regression-2d-flower-objects.rds"
))
objects <- upgrade_cv_mse(objects)
z_external <- objects$data$z_external
response_observed <- objects$data$response_observed
conditional_mean <- objects$data$conditional_mean
fold_id <- objects$data$fold_id
lambda_grid <- objects$tuning$lambda_grid
bandwidth_grid <- objects$tuning$bandwidth_grid
if (identical(Sys.getenv("PME_SMOKE_TEST"), "1") && exists("bandwidth_grid")) {
  bandwidth_grid <- c(0.2, 0.4, 0.8)
}
methods <- c("kernel", "local_linear")
surface_grid <- fibonaccisphere(if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 100L else 1800L)

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
    mean_raw_mse = mean(group$raw_rmse^2),
    se_raw_mse = sd(group$raw_rmse^2) / sqrt(nrow(group)),
    mean_oracle_mse = mean(group$oracle_rmse^2),
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
  x[which.min(x$mean_raw_mse), , drop = FALSE]
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

}
if (identical(Sys.getenv("PME_SMOKE_TEST"), "1") ||
    !file.exists(file.path(result_dir, "fixed-bandwidth-local-linear-objects.rds"))) {
  if (!all(file.exists(file.path(cache_dir, sprintf(
    "pme-lambda-%02d-fold-%d.rds", rep(seq_along(readRDS(input_file)$tuning$lambda_grid), each = 5L),
    rep(seq_len(5L), length(readRDS(input_file)$tuning$lambda_grid))))))) reproduce_inputs(notebook_dir)
  reproduce_fixed_bandwidth()
}
if (!file.exists(file.path(result_dir, "bandwidth-sensitivity-all-combinations.csv"))) {
  if (!all(file.exists(file.path(cache_dir, sprintf(
    "pme-lambda-%02d-fold-%d.rds", rep(seq_along(readRDS(input_file)$tuning$lambda_grid), each = 5L),
    rep(seq_len(5L), length(readRDS(input_file)$tuning$lambda_grid))))))) reproduce_inputs(notebook_dir)
  reproduce_sensitivity()
}
objects <- readRDS(file.path(
  result_dir, "fixed-bandwidth-local-linear-objects.rds"
))
objects <- upgrade_cv_mse(objects)
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
if (identical(Sys.getenv("PME_SMOKE_TEST"), "1") && exists("bandwidth_grid")) {
  bandwidth_grid <- c(0.2, 0.4, 0.8)
}
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
  best <- values[which.min(values$mean_raw_mse), ]
  y_range <- range(
    values$mean_raw_mse + c(-1, 1) * values$se_raw_mse,
    values$mean_oracle_mse
  )
  plot(values$lambda, values$mean_raw_mse,
       type = "b", log = "x", pch = 16, lwd = 2.2,
       col = "#D55E00", ylim = y_range,
       xlab = expression(lambda), ylab = "Five-fold MSE",
       main = sprintf("%s, h = %.2f", labels[method], objects$bandwidth))
  arrows(values$lambda,
         values$mean_raw_mse - values$se_raw_mse,
         values$lambda,
         values$mean_raw_mse + values$se_raw_mse,
         angle = 90, code = 3, length = 0.035, col = "#D55E00")
  lines(values$lambda, values$mean_oracle_mse,
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
