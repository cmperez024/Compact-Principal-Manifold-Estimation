## Regenerate the comparison figure from saved objects without refitting PME.

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
local_library <- file.path(repo_root, "packages")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
if (!requireNamespace("pracma", quietly = TRUE)) stop("pracma is required.")
source(file.path(repo_root, "CompactPME", "R", "krt.R"))
source(file.path(repo_root, "CompactPME", "R", "spline1d_periodic.R"))
bernoulli <- pracma::bernoulli
ones <- pracma::ones
eye <- pracma::eye
result_dir <- file.path(manuscript_experiments, experiment_name, "results")
objects <- readRDS(file.path(
  result_dir, "fixed-design-extrinsic-regression-objects.rds"
))
response_observed <- objects$data$response_observed
true_flower <- objects$display$true_mean_grid
final <- objects$final
methods <- c("kernel", "local_linear")
colors <- c(kernel = "#0072B2", local_linear = "#009E73")
titles <- c(kernel = "Kernel-selected PME",
            local_linear = "Local-linear-selected PME")

png(file.path(result_dir, "01-fixed-design-kernel-vs-local-linear.png"),
    width = 3200, height = 1200, res = 210)
par(mfrow = c(1, 3), mar = c(4.3, 4.3, 3.4, 1), mgp = c(2.5, 0.75, 0))
plot(response_observed, asp = 1, pch = 16, cex = 0.48,
     col = adjustcolor("grey35", 0.32), xlab = expression(X[1]),
     ylab = expression(X[2]), main = "Fixed-design response-noise data")
lines(true_flower, col = "black", lwd = 2.5)
legend("topright", c("observed response", "conditional mean"),
       pch = c(16, NA), lty = c(NA, 1), col = c("grey45", "black"),
       lwd = c(NA, 2.5), bty = "n")

for (method in methods) {
  choice <- final[[method]]$choice
  pme_curve <- final[[method]]$pme$fit(seq(0, 1, length.out = 2001L))
  plot(response_observed, asp = 1, pch = 16, cex = 0.42,
       col = adjustcolor("grey45", 0.22), xlab = expression(X[1]),
       ylab = expression(X[2]),
       main = sprintf("%s\nlambda=%.0e",
                      titles[method], choice$lambda))
  lines(true_flower, col = "black", lwd = 2.5, lty = 2)
  lines(pme_curve, col = colors[method], lwd = 3)
  legend("topright", c("true flower", "selected PME manifold"),
         col = c("black", colors[method]), lty = c(2, 1),
         lwd = c(2.5, 3), bty = "n")
}
dev.off()
