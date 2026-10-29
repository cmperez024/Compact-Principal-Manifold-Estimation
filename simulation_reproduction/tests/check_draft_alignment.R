# Lightweight regression checks; no simulations or publication outputs are run.
script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script[1L])), "..", ".."))
renderer <- file.path(root, "simulation_reproduction", "01_main_text",
  "04_downstream_tuning", "code_and_results", "2d_flower",
  "render-local-linear-composite-figures.R")
expressions <- parse(renderer)
for (e in expressions) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      identical(e[[2]], as.name("upgrade_cv_mse"))) eval(e)
}
folds <- data.frame(method = "local_linear", lambda = c(1, 1, 2, 2),
  fold = c(1, 2, 1, 2), validation_rmse_raw = c(0.1, 0.9, 0.55, 0.55),
  validation_rmse_oracle = c(0.1, 0.9, 0.55, 0.55))
scores <- data.frame(method = "local_linear", lambda = c(1, 2),
                     mean_raw_rmse = c(0.5, 0.55), selected = c(TRUE, FALSE))
x <- upgrade_cv_mse(list(tuning = list(fold_scores = folds,
  aggregate_scores = scores, selected = list()), final = list()))
stopifnot(x$tuning$selected$local_linear$lambda == 2,
          isTRUE(all.equal(x$tuning$aggregate_scores$mean_raw_mse,
                           c(0.41, 0.3025))),
          identical(x$tuning$selection_loss, "mean_fold_mse"))
folds$validation_rmse_raw[1] <- 0.2
legacy <- list(tuning = list(fold_scores = folds, aggregate_scores = scores,
  selected = list()), final = list(local_linear = list(choice = list(lambda = 1))))
stopifnot(inherits(try(upgrade_cv_mse(legacy), silent = TRUE), "try-error"))
legacy$tuning$fold_scores$lambda <- c(1e-8, 1e-8, 1e-9, 1e-9)
legacy$tuning$aggregate_scores$lambda <- c(1e-8, 1e-9)
legacy$final$local_linear$choice$lambda <- 1e-8
# Distinct tiny lambdas must not be accepted by all.equal's default tolerance.
stopifnot(inherits(try(upgrade_cv_mse(legacy), silent = TRUE), "try-error"))
cat("Draft-alignment regression checks passed.\n")
