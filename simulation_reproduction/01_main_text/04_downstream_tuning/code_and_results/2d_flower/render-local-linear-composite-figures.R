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

## Paste-ready local-linear criterion + truth + selected-PME figures for the
## current two-dimensional flower experiment and the preceding one-dimensional
## flower experiment.

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
if (!requireNamespace("CompactPME", quietly = TRUE)) {
  stop("Install the CompactPME package before running this script.")
}
qm <- CompactPME:::qm
expr <- CompactPME:::expr
krt <- CompactPME:::krt
krst <- CompactPME:::krst
Rk <- CompactPME:::Rk
spline2d <- CompactPME::spline2d
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
draw_criterion <- function(lambda, raw, selected_lambda,
                           selected_raw, main, bandwidth_text) {
  y_limits <- range(raw)
  plot(
    lambda, raw, type = "b", log = "x", pch = 16,
    lwd = 2.4, col = "#D55E00", ylim = y_limits,
    xlab = expression(lambda), ylab = "Five-fold out-of-fold MSE",
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

# This figure owns both experiments and can regenerate missing inputs.
# Generator bodies are included below; no external simulation is executed.
reproduce_1d <- function(target_dir) {
options(stringsAsFactors = FALSE)
set.seed(20260806)

notebook_dir <- target_dir
result_dir <- file.path(notebook_dir, "results")
cache_dir <- file.path(notebook_dir, "cache")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

local_library <- file.path(repo_root, "packages")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
if (!requireNamespace("CompactPME", quietly = TRUE)) {
  stop("Install the CompactPME package before running this script.")
}
spline1d_periodic <- CompactPME::spline1d_periodic
project_circle <- CompactPME:::project_circle
project_grid <- CompactPME:::project_grid

flower_map <- function(u, petals = 5L, amplitude = 0.30) {
  angle <- 2 * pi * u
  radius <- 1 + amplitude * cos(petals * angle)
  cbind(radius * cos(angle), radius * sin(angle))
}

phase_map <- function(t) (0.5 + 0.5 * sin(pi * t)) %% 1

fit_pme_s1 <- function(X, lambda, projection_grid_size = if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 100L else 1200L,
                       min_iterations = 3L, max_iterations = if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 3L else 100L,
                       relative_tolerance = 1e-4) {
  projection_index <- project_circle(X)
  fit <- spline1d_periodic(X, projection_index, lambda)
  projection_index <- project_grid(X, fit, gridSize = projection_grid_size)
  fidelity <- mean(rowSums((X - fit(projection_index))^2))
  history <- data.frame(iteration = 1L, fidelity = fidelity,
                        relative_change = NA_real_)
  converged <- FALSE
  for (iteration in 2:max_iterations) {
    next_fit <- spline1d_periodic(X, projection_index, lambda)
    next_index <- project_grid(X, next_fit, gridSize = projection_grid_size)
    next_fidelity <- mean(rowSums((X - next_fit(next_index))^2))
    relative_change <- abs(next_fidelity - fidelity) / max(fidelity, 1e-14)
    history <- rbind(history, data.frame(
      iteration = iteration, fidelity = next_fidelity,
      relative_change = relative_change
    ))
    fit <- next_fit
    projection_index <- next_index
    fidelity <- next_fidelity
    if (iteration >= min_iterations && relative_change <= relative_tolerance) {
      converged <- TRUE
      break
    }
  }
  list(fit = fit, projection_index = projection_index,
       denoised = fit(projection_index), history = history,
       converged = converged, lambda = lambda)
}

ambient_kernel <- function(t_train, y_train, t_new, bandwidth,
                           leave_one_out = FALSE) {
  weights <- exp(-0.5 * (outer(t_new, t_train, "-") / bandwidth)^2)
  if (leave_one_out) {
    if (length(t_train) != length(t_new) ||
        !isTRUE(all.equal(t_train, t_new))) {
      stop("LOO kernel prediction requires identical training and new t.")
    }
    diag(weights) <- 0
  }
  weight_sums <- rowSums(weights)
  if (any(weight_sums <= .Machine$double.eps)) {
    stop("Bandwidth produced effectively zero total kernel weight.")
  }
  weights %*% y_train / weight_sums
}

ambient_local_linear <- function(t_train, y_train, t_new, bandwidth,
                                 leave_one_out = FALSE) {
  if (leave_one_out &&
      (length(t_train) != length(t_new) ||
       !isTRUE(all.equal(t_train, t_new)))) {
    stop("LOO local-linear prediction requires identical training and new t.")
  }
  prediction <- matrix(NA_real_, nrow = length(t_new), ncol = ncol(y_train))
  for (j in seq_along(t_new)) {
    z <- t_train - t_new[j]
    w <- exp(-0.5 * (z / bandwidth)^2)
    if (leave_one_out) w[j] <- 0
    s0 <- sum(w)
    if (s0 <= .Machine$double.eps) {
      stop("Bandwidth produced effectively zero total local-linear weight.")
    }
    s1 <- sum(w * z)
    s2 <- sum(w * z^2)
    denominator <- s0 * s2 - s1^2
    if (denominator <= 1e-12 * max(1, s0 * s2)) {
      prediction[j, ] <- colSums(w * y_train) / s0
    } else {
      t0 <- colSums(w * y_train)
      t1 <- colSums((w * z) * y_train)
      prediction[j, ] <- (s2 * t0 - s1 * t1) / denominator
    }
  }
  prediction
}

extrinsic_predict <- function(method, t_train, y_train, t_new, bandwidth,
                              manifold_fit, projection_grid_size = if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 100L else 1200L,
                              leave_one_out = FALSE) {
  ambient <- switch(method,
    kernel = ambient_kernel(
      t_train, y_train, t_new, bandwidth, leave_one_out
    ),
    local_linear = ambient_local_linear(
      t_train, y_train, t_new, bandwidth, leave_one_out
    )
  )
  index <- project_grid(ambient, manifold_fit$fit,
                        gridSize = projection_grid_size)
  list(ambient = ambient, index = index,
       extrinsic = manifold_fit$fit(index))
}

rmse_2d <- function(observed, predicted) {
  sqrt(mean(rowSums((observed - predicted)^2)))
}

select_bandwidth <- function(method, t_train, denoised_train,
                             bandwidth_grid, manifold_fit) {
  error <- vapply(bandwidth_grid, function(bandwidth) {
    prediction <- extrinsic_predict(
      method, t_train, denoised_train, t_train, bandwidth,
      manifold_fit, leave_one_out = TRUE
    )$extrinsic
    rmse_2d(denoised_train, prediction)
  }, numeric(1))
  selected <- which.min(error)
  list(
    bandwidth = bandwidth_grid[selected],
    rmse = error[selected],
    errors = error
  )
}

## Generate the covariates once. From here onward they are conditioned upon.
n <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 40L else 600L
noise_sd <- as.numeric(Sys.getenv("PME_RESPONSE_NOISE_SD", unset = "0.04"))
if (!is.finite(noise_sd) || noise_sd <= 0) stop("Invalid response-noise SD.")
t_fixed <- sort(runif(n, -1, 1))
u_deterministic <- phase_map(t_fixed)
conditional_mean <- flower_map(u_deterministic)
response_error <- matrix(rnorm(2L * n, sd = noise_sd), ncol = 2L)
response_observed <- conditional_mean + response_error

## Generate the validation partition once and hold it fixed as well.
k_folds <- 5L
fold_id <- sample(rep(seq_len(k_folds), length.out = n))
lambda_grid <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) c(1e-6, 1e-4) else 10^seq(-14, -4, by = 1)
bandwidth_grid <- c(0.004, 0.006, 0.008, 0.01, 0.015, 0.022, 0.03, 0.045, 0.065, 0.09,
                    0.13, 0.18, 0.26, 0.38)
if (identical(Sys.getenv("PME_SMOKE_TEST"), "1") && exists("bandwidth_grid")) {
  bandwidth_grid <- c(0.2, 0.4, 0.8)
}
methods <- c("kernel", "local_linear")
prediction_store <- array(
  NA_real_, dim = c(n, 2L, length(lambda_grid), length(methods)),
  dimnames = list(NULL, c("X1", "X2"), paste(lambda_grid),
                  methods)
)
fit_rows <- vector("list", length(lambda_grid) * k_folds)
score_rows <- vector("list", length(lambda_grid) * length(methods) * k_folds)
fit_row <- 0L
score_row <- 0L

for (li in seq_along(lambda_grid)) {
  lambda <- lambda_grid[li]
  for (fold in seq_len(k_folds)) {
    validation_id <- which(fold_id == fold)
    training_id <- which(fold_id != fold)
    cache_path <- file.path(
      cache_dir, sprintf("pme-lambda-%02d-fold-%d.rds", li, fold)
    )
    if (file.exists(cache_path)) {
      pme <- readRDS(cache_path)
      if (!isTRUE(pme$converged)) {
        pme <- fit_pme_s1(
          response_observed[training_id, , drop = FALSE], lambda
        )
        saveRDS(pme, cache_path)
      }
    } else {
      pme <- fit_pme_s1(response_observed[training_id, , drop = FALSE], lambda)
      saveRDS(pme, cache_path)
    }
    fit_row <- fit_row + 1L
    fit_rows[[fit_row]] <- data.frame(
      lambda = lambda, fold = fold, iterations = nrow(pme$history),
      converged = pme$converged
    )
    for (mi in seq_along(methods)) {
      method <- methods[mi]
      bandwidth_selection <- select_bandwidth(
        method, t_fixed[training_id], pme$denoised, bandwidth_grid, pme
      )
      prediction <- extrinsic_predict(
        method, t_fixed[training_id], pme$denoised,
        t_fixed[validation_id], bandwidth_selection$bandwidth, pme
      )$extrinsic
      prediction_store[validation_id, , li, mi] <- prediction
      score_row <- score_row + 1L
      score_rows[[score_row]] <- data.frame(
        method = method,
        lambda = lambda,
        fold = fold,
        bandwidth = bandwidth_selection$bandwidth,
        training_loo_rmse = bandwidth_selection$rmse,
        validation_rmse_raw = rmse_2d(
          response_observed[validation_id, , drop = FALSE], prediction
        ),
        validation_rmse_oracle = rmse_2d(
          conditional_mean[validation_id, , drop = FALSE], prediction
        )
      )
    }
    cat(sprintf("lambda=%g fold=%d iterations=%d converged=%s\n",
                lambda, fold, nrow(pme$history), pme$converged))
  }
}

fold_scores <- do.call(rbind, score_rows)
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
  fold_scores,
  interaction(fold_scores$method, fold_scores$lambda, drop = TRUE)
)
aggregate_scores <- do.call(rbind, lapply(groups, function(group) {
  cbind(
    data.frame(method = group$method[1], lambda = group$lambda[1]),
    aggregate_one(group)
  )
}))
rownames(aggregate_scores) <- NULL
aggregate_scores <- aggregate_scores[
  order(aggregate_scores$method, aggregate_scores$lambda),
]
aggregate_scores$selected <- FALSE
selected <- vector("list", length(methods))
names(selected) <- methods
for (method in methods) {
  eligible <- which(aggregate_scores$method == method)
  best <- eligible[which.min(aggregate_scores$mean_raw_mse[eligible])]
  aggregate_scores$selected[best] <- TRUE
  selected[[method]] <- aggregate_scores[best, ]
}

## Refit each selected PME on all responses and produce smooth display curves.
t_grid <- seq(-1, 1, length.out = 501L)
true_mean_grid <- flower_map(phase_map(t_grid))
final <- vector("list", length(methods))
names(final) <- methods
full_fit_by_lambda <- list()
for (method in methods) {
  choice <- selected[[method]]
  lambda <- choice$lambda
  cache_key <- format(lambda, scientific = TRUE, digits = 16)
  if (is.null(full_fit_by_lambda[[cache_key]])) {
    full_fit_by_lambda[[cache_key]] <- fit_pme_s1(response_observed, lambda)
  }
  pme <- full_fit_by_lambda[[cache_key]]
  bandwidth_selection <- select_bandwidth(
    method, t_fixed, pme$denoised, bandwidth_grid, pme
  )
  prediction <- extrinsic_predict(
    method, t_fixed, pme$denoised, t_grid,
    bandwidth_selection$bandwidth, pme
  )
  final[[method]] <- list(
    choice = choice,
    bandwidth = bandwidth_selection$bandwidth,
    bandwidth_selection = bandwidth_selection,
    pme = pme,
    t_grid = t_grid,
    prediction = prediction,
    oof_prediction = prediction_store[, ,
      which(lambda_grid == lambda), which(methods == method)]
  )
}

write.csv(fold_scores, file.path(result_dir, "fold-performance.csv"),
          row.names = FALSE)
write.csv(aggregate_scores, file.path(result_dir, "tuning-performance.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, fit_rows), file.path(result_dir, "pme-fit-status.csv"),
          row.names = FALSE)
saveRDS(list(
  model = list(n = n, noise_sd = noise_sd,
               phase_map = "u(t) = (1 + sin(pi*t))/2 modulo 1"),
  data = list(t_fixed = t_fixed, u_deterministic = u_deterministic,
              conditional_mean = conditional_mean,
              response_error = response_error,
              response_observed = response_observed, fold_id = fold_id),
  tuning = list(lambda_grid = lambda_grid, bandwidth_grid = bandwidth_grid,
                fold_scores = fold_scores,
                aggregate_scores = aggregate_scores, selected = selected),
  predictions = prediction_store,
  final = final,
  display = list(t_grid = t_grid, true_mean_grid = true_mean_grid)
), file.path(result_dir, "fixed-design-extrinsic-regression-objects.rds"))

png(file.path(result_dir, "01-fixed-design-kernel-vs-local-linear.png"),
    width = 3200, height = 1200, res = 210)
par(mfrow = c(1, 3), mar = c(4.3, 4.3, 3.4, 1), mgp = c(2.5, 0.75, 0))
plot(response_observed, asp = 1, pch = 16, cex = 0.48,
     col = adjustcolor("grey35", 0.32), xlab = expression(X[1]),
     ylab = expression(X[2]), main = "Fixed-design response-noise data")
lines(true_mean_grid, col = "black", lwd = 2.5)
legend("topright", c("observed response", "conditional mean"),
       pch = c(16, NA), lty = c(NA, 1), col = c("grey45", "black"),
       lwd = c(NA, 2.5), bty = "n")

colors <- c(kernel = "#0072B2", local_linear = "#009E73")
titles <- c(kernel = "Kernel-selected PME", local_linear = "Local-linear-selected PME")
for (method in methods) {
  choice <- selected[[method]]
  pme_curve <- final[[method]]$pme$fit(seq(0, 1, length.out = 2001L))
  plot(response_observed, asp = 1, pch = 16, cex = 0.42,
       col = adjustcolor("grey45", 0.22), xlab = expression(X[1]),
       ylab = expression(X[2]),
       main = sprintf("%s\nlambda=%.0e",
                      titles[method], choice$lambda))
  lines(true_mean_grid, col = "black", lwd = 2.5, lty = 2)
  lines(pme_curve, col = colors[method], lwd = 3)
  legend("topright", c("true flower", "selected PME manifold"),
         col = c("black", colors[method]), lty = c(2, 1),
         lwd = c(2.5, 3), bty = "n")
}
dev.off()

for (method in methods) {
  choice <- selected[[method]]
  cat(sprintf("%s selected lambda=%.0e h=%.3f OOF noisy MSE=%.5f oracle MSE=%.5f\n",
              method, choice$lambda, final[[method]]$bandwidth,
              choice$mean_raw_mse, choice$mean_oracle_mse))
}

}
reproduce_2d <- function(target_dir, variant = "uniform-noninjective", n_sample = 800L) {
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
      "%s OOF predictions\nraw MSE = %.4f",
      method_labels[method], final[[method]]$choice$mean_raw_mse
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
    "%s: lambda=%.3e, mean OOF raw MSE=%.5f, oracle MSE=%.5f, full h=%.2f\n",
    method_labels[method], choice$lambda, choice$mean_raw_mse,
    choice$mean_oracle_mse, final[[method]]$bandwidth
  ))
}

}
figure_inputs <- Sys.getenv("PME_RESULTS_ROOT", unset = file.path(manuscript_section, "figure_5_2_inputs"))
dir.create(figure_inputs, recursive = TRUE, showWarnings = FALSE)
notebook_1d <- file.path(figure_inputs, "1d")
notebook_2d <- file.path(figure_inputs, "2d")
if (!file.exists(file.path(notebook_1d, "results",
                          "fixed-design-extrinsic-regression-objects.rds"))) {
  reproduce_1d(notebook_1d)
}
if (!file.exists(file.path(notebook_2d, "results-uniform-n800",
                          "extrinsic-regression-2d-flower-objects.rds"))) {
  reproduce_2d(notebook_2d)
}

## Two-dimensional flower --------------------------------------------------
result_2d <- file.path(notebook_2d, "results-uniform-n800")
objects_2d <- readRDS(file.path(
  result_2d, "extrinsic-regression-2d-flower-objects.rds"
))
objects_2d <- upgrade_cv_mse(objects_2d)
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
result_1d <- file.path(notebook_1d, "results")
drop_1d <- file.path(result_1d, "drop-in")
objects_1d <- readRDS(file.path(
  result_1d, "fixed-design-extrinsic-regression-objects.rds"
))
objects_1d <- upgrade_cv_mse(objects_1d)
profile_1d <- objects_1d$tuning$aggregate_scores
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
    profile_1d$lambda, profile_1d$mean_raw_mse,
    selected_1d$lambda, selected_1d$mean_raw_mse, title,
    sprintf("mean fold-selected bandwidth = %.3f",
            selected_1d$mean_bandwidth)
  )
}

draw_criterion_2d <- function(title) {
  draw_criterion(
    profile_2d$lambda, profile_2d$mean_raw_mse,
    selected_2d$lambda, selected_2d$mean_raw_mse, title,
    sprintf("mean fold-selected bandwidth = %.3f", selected_2d$mean_bandwidth)
  )
}

## Combined 2 x 3 publication figure. The 1D fit spans the upper-right cells.
## Write only the figure referenced by the manuscript; component renderings are
## deliberately not emitted into the reproduction tree.
output_dir <- Sys.getenv("PME_OUTPUT_DIR", unset = manuscript_section)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
paper_output_combined <- file.path(output_dir, "figure_5_2_downstream_tuning.png")
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
