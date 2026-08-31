## Publication-ready figures for the fixed-design SD 0.04 experiment.

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
notebook_dir <- file.path(
  manuscript_experiments,
  "1d_sd04"
)
result_dir <- file.path(notebook_dir, "results")
drop_dir <- file.path(result_dir, "drop-in")
dir.create(drop_dir, recursive = TRUE, showWarnings = FALSE)

local_library <- file.path(repo_root, "packages")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
if (!requireNamespace("pracma", quietly = TRUE)) stop("pracma is required.")
source(file.path(repo_root, "CompactPME", "R", "krt.R"))
source(file.path(repo_root, "CompactPME", "R", "spline1d_periodic.R"))
bernoulli <- pracma::bernoulli
ones <- pracma::ones
eye <- pracma::eye

objects <- readRDS(file.path(
  result_dir, "fixed-design-extrinsic-regression-objects.rds"
))
scores <- objects$tuning$aggregate_scores
methods <- c("kernel", "local_linear")
method_titles <- c(kernel = "Extrinsic kernel regression",
                   local_linear = "Extrinsic local-linear regression")

## Each row already uses the fold-specific bandwidths selected by LOO on the
## corresponding training fold. The outer validation responses select lambda.
profile <- scores

png(file.path(drop_dir, "kfold_performance.png"),
    width = 2600, height = 1200, res = 220)
par(mfrow = c(1, 2), mar = c(4.5, 4.7, 3.3, 1.0),
    mgp = c(2.7, 0.8, 0))
all_y <- range(c(profile$mean_raw_rmse, profile$mean_oracle_rmse))
for (method in methods) {
  z <- profile[profile$method == method, ]
  selected <- objects$tuning$selected[[method]]
  plot(z$lambda, z$mean_raw_rmse, type = "b", log = "x", pch = 16,
       lwd = 2.4, col = "#D55E00", ylim = all_y,
       xlab = expression(lambda), ylab = "Five-fold out-of-fold RMSE",
       main = method_titles[method])
  lines(z$lambda, z$mean_oracle_rmse, type = "b", pch = 1,
        lwd = 2.4, col = "#0072B2", lty = 2)
  abline(v = selected$lambda, col = "grey35", lwd = 2, lty = 3)
  points(selected$lambda, selected$mean_raw_rmse, pch = 21, cex = 1.45,
         bg = "#D55E00", col = "black", lwd = 1.2)
  legend("topleft",
         c("untouched noisy response", "conditional mean (oracle)",
           sprintf("selected: lambda=%.0e", selected$lambda)),
         col = c("#D55E00", "#0072B2", "grey35"),
         lty = c(1, 2, 3), pch = c(16, 1, NA), lwd = 2.2,
         bty = "n", cex = 0.84)
  mtext(sprintf("mean fold-selected bandwidth = %.3f",
                selected$mean_bandwidth), side = 3, line = 0.25,
        cex = 0.72, col = "grey35")
}
dev.off()

response <- objects$data$response_observed
truth <- objects$display$true_mean_grid
selected_lambda <- objects$tuning$selected$local_linear$lambda
selected_pme <- objects$final$local_linear$pme$fit(
  seq(0, 1, length.out = 2001L)
)

png(file.path(drop_dir, "kfold_fits.png"),
    width = 2600, height = 1200, res = 220)
par(mfrow = c(1, 2), mar = c(4.5, 4.6, 3.3, 1.0),
    mgp = c(2.7, 0.8, 0))
plot(response, asp = 1, pch = 16, cex = 0.48,
     col = adjustcolor("grey40", alpha.f = 0.32),
     xlab = expression(X[1]), ylab = expression(X[2]),
     main = "Fixed-design responses")
lines(truth, col = "black", lwd = 2.7)
legend("topright", c("observed response", "true conditional-mean curve"),
       pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2.7),
       col = c("grey45", "black"), bty = "n", cex = 0.85)

plot(response, asp = 1, pch = 16, cex = 0.42,
     col = adjustcolor("grey50", alpha.f = 0.22),
     xlab = expression(X[1]), ylab = expression(X[2]),
     main = sprintf("Cross-validation-selected PME: lambda = %.0e",
                    selected_lambda))
lines(truth, col = "black", lwd = 2.5, lty = 2)
lines(selected_pme, col = "#0072B2", lwd = 3)
legend("topright", c("true flower", "selected PME manifold"),
       col = c("black", "#0072B2"), lty = c(2, 1),
       lwd = c(2.5, 3), bty = "n", cex = 0.85)
dev.off()

write.csv(profile, file.path(drop_dir, "lambda-performance-profile.csv"),
          row.names = FALSE)
