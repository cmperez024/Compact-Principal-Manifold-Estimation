## Fixed-design response-noise experiment for Lin et al. (2017): compare
## extrinsic kernel and extrinsic local-linear regression after PME.

options(stringsAsFactors = FALSE)
set.seed(20260806)

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
manuscript_experiments <- file.path(
  repo_root, "simulation_reproduction", "01_main_text",
  "04_downstream_tuning", "code_and_results"
)
experiment_name <- Sys.getenv(
  "PME_FIXED_DESIGN_EXPERIMENT",
  unset = "1d_base"
)
notebook_dir <- file.path(manuscript_experiments, experiment_name)
result_dir <- file.path(notebook_dir, "results")
cache_dir <- file.path(notebook_dir, "cache")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

local_library <- file.path(repo_root, "packages")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
if (!requireNamespace("pracma", quietly = TRUE)) stop("pracma is required.")
source(file.path(repo_root, "CompactPME", "R", "krt.R"))
source(file.path(repo_root, "CompactPME", "R", "spline1d_periodic.R"))
source(file.path(repo_root, "CompactPME", "R", "project_circle.R"))
bernoulli <- pracma::bernoulli
ones <- pracma::ones
eye <- pracma::eye

flower_map <- function(u, petals = 5L, amplitude = 0.30) {
  angle <- 2 * pi * u
  radius <- 1 + amplitude * cos(petals * angle)
  cbind(radius * cos(angle), radius * sin(angle))
}

phase_map <- function(t) (0.5 + 0.5 * sin(pi * t)) %% 1

fit_pme_s1 <- function(X, lambda, projection_grid_size = 1200L,
                       min_iterations = 3L, max_iterations = 100L,
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
                              manifold_fit, projection_grid_size = 1200L,
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
n <- 600L
noise_sd <- as.numeric(Sys.getenv("PME_RESPONSE_NOISE_SD", unset = "0.08"))
if (!is.finite(noise_sd) || noise_sd <= 0) stop("Invalid response-noise SD.")
t_fixed <- sort(runif(n, -1, 1))
u_deterministic <- phase_map(t_fixed)
conditional_mean <- flower_map(u_deterministic)
response_error <- matrix(rnorm(2L * n, sd = noise_sd), ncol = 2L)
response_observed <- conditional_mean + response_error

## Generate the validation partition once and hold it fixed as well.
k_folds <- 5L
fold_id <- sample(rep(seq_len(k_folds), length.out = n))
lambda_grid <- 10^seq(-14, -4, by = 1)
bandwidth_grid <- c(0.004, 0.006, 0.008, 0.01, 0.015, 0.022, 0.03, 0.045, 0.065, 0.09,
                    0.13, 0.18, 0.26, 0.38)
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
  best <- eligible[which.min(aggregate_scores$mean_raw_rmse[eligible])]
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
  cat(sprintf("%s selected lambda=%.0e h=%.3f OOF noisy RMSE=%.5f oracle RMSE=%.5f\n",
              method, choice$lambda, final[[method]]$bandwidth,
              choice$mean_raw_rmse, choice$mean_oracle_rmse))
}
