# Reproduce Figure 4.1: empirical consistency of the periodic PME.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1L])))

if (!requireNamespace("CompactPME", quietly = TRUE)) {
  stop("Install CompactPME before running this script.")
}

sample_sizes <- c(100L, 250L, 500L, 1000L)
repetitions <- 100L
noise_sd <- 0.04
lambda <- 1e-8
evaluation_grid <- seq(0, 1, length.out = 500L)
set.seed(1)

flower_curve <- function(t, noise = 0) {
  radius <- 1 + 0.3 * sin(5 * 2 * pi * t)
  curve <- cbind(radius * cos(2 * pi * t), radius * sin(2 * pi * t))
  curve + matrix(rnorm(2L * length(t), sd = noise), ncol = 2L)
}

fit_once <- function(n) {
  observations <- flower_curve(sort(runif(n)), noise_sd)
  fit <- CompactPME::pme_1d_periodic(
    observations,
    lambdas = lambda,
    optimize_lambda = FALSE,
    parallel = FALSE,
    approx = FALSE
  )
  fit$spline_list[[1L]](evaluation_grid)
}

fits <- lapply(sample_sizes, function(n) {
  replicate(repetitions, fit_once(n), simplify = FALSE)
})
truth <- flower_curve(evaluation_grid)
limit <- max(abs(unlist(c(list(truth), fits))))

default_output <- file.path(
  dirname(script_dir), "PME_monte_carlo_star2.png"
)
output_file <- Sys.getenv("PME_FIGURE_OUTPUT", unset = default_output)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

png(output_file, width = 12, height = 4, units = "in", res = 600)
layout(
  matrix(c(1:4, rep(5, 4)), nrow = 2, byrow = TRUE),
  heights = c(1, 0.14)
)
for (i in seq_along(sample_sizes)) {
  plot(
    truth,
    type = "n",
    asp = 1,
    xlim = c(-limit, limit),
    ylim = c(-limit, limit),
    xlab = "x",
    ylab = "y",
    main = paste("N =", sample_sizes[i])
  )
  for (curve in fits[[i]]) {
    lines(curve, col = grDevices::adjustcolor("blue", alpha.f = 0.3))
  }
  lines(truth, col = "red", lwd = 2)
}
par(mar = rep(0, 4))
plot.new()
legend(
  "center",
  legend = c("Sample Curves", "True Curve"),
  col = c("blue", "red"),
  lty = 1,
  horiz = TRUE,
  bty = "o"
)
dev.off()

cat("Wrote", normalizePath(output_file, winslash = "/", mustWork = TRUE), "\n")
cat("lambda =", format(lambda, scientific = TRUE), "\n")
