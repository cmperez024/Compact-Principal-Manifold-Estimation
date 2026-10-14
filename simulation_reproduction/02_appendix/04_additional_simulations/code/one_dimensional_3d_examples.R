# Resolve paths from this file, including when invoked outside the project.
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_file)) stop("Run this script with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_file[1L]),
                                     winslash = "/", mustWork = TRUE))
figure_dir <- Sys.getenv("PME_OUTPUT_DIR", unset = file.path(script_dir, ".."))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- normalizePath(figure_dir, mustWork = TRUE)
repo_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."),
                           mustWork = TRUE)
local_library <- file.path(repo_root, "packages")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
setwd(figure_dir)
set.seed(1)
for (folder in c("results_1d", "results_sph", "job_figs_1d",
                  "panel_flower", "panel_cashew", "sim_results")) {
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
}

# Reproduce only the two appendix 1D-in-3D figures.
# All data generators are included here; fitting comes from CompactPME.
required <- c("CompactPME", "plot3D", "wacolors")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " first.")
}
library(plot3D)
library(wacolors)
smoke <- identical(Sys.getenv("PME_SMOKE_TEST"), "1")
n_sample <- if (smoke) 40L else 500L
fit_iterations <- if (smoke) 2L else 10L
lambdas <- 10^c(-12, -6, 0)

rotation3 = function(u, angle=NULL) {
  n = length(u)
  if (n==3) {
    unorm = sqrt(sum(u^2))
    if (unorm==0) {
      R = diag(3)
    } else {
      if (is.null(angle)) angle = unorm
      u = u/unorm
      c = cos(angle)
      s = sin(angle)
      R = rbind(c(c+u[1]^2*(1-c), u[1]*u[2]*(1-c)-u[3]*s, u[1]*u[3]*(1-c)+u[2]*s),
                c(u[2]*u[1]*(1-c)+u[3]*s, c+u[2]^2*(1-c), u[2]*u[3]*(1-c)-u[1]*s),
                c(u[3]*u[1]*(1-c)-u[2]*s, u[3]*u[2]*(1-c)+u[1]*s, c+u[3]^2*(1-c)))
    }
    
    return(R)
  } else {
    stop('u must be a 3-element vector')
  }
}

sine1d3D <- function(n, rotation=0, period=2*pi, noise=0, tval){
  t0 <- if(missing(tval)) sort(runif(n)) else tval
  X <- t(cbind(t0, sin(period*t0),cos(period*t0)))
  X <- rotation3(c(1,1,0.2),rotation)%*% X
  X <- t(X) + rnorm(3*n, sd=noise)
  list(X=X,t=t0)
}

helix1d3D <- function(n, period = 4*pi, xyscale = 1, zscale=1, noise=0, rotation=0){
  
  t0 <- sort(runif(n))
  
  list(X=cbind(xyscale*cos(period*t0), xyscale*sin(period*t0), zscale*t0) %*% t( rotation3(c(1,1,0.2),rotation))+ rnorm(3*n, sd=noise),
       t = t0)
}

flower1d3D <- function(n, petals, noise, t){
  t0 <- if(missing(t)) sort(runif(n)) else t
  
  r <- (1 + 0.3 * sin(petals *2*pi* t0))
  return(list(t=t0, X = cbind(r * cos(2*pi*t0),  r * sin(2*pi*t0), sin(3*pi*t0))+rnorm(3*n, sd=noise)))
}

moon1d3D <- function(N, inner_radius=1, outer_radius=2, noise=0)
{
  # Parameters (example)
  r1 <- inner_radius    # inner radius
  r2 <- outer_radius      # outer radius (must be > r1)
  stopifnot(r2 > r1)
  #n <- N/4 # points per segment (increase for smoother curve)
  nmain <- round(4*N/5)
  ncap <- round(N/5)
  
  
  r_cap <- (r2 - r1) / 2
  k <- (r1 + r2) / 2    # midpoint y for cap centers
  
  # 1) Outer left semicircle: top -> bottom (theta from pi/2 to 3pi/2)
  theta_outer <- seq(pi/2, 3*pi/2, length.out = nmain)
  x_outer <- r2 * cos(theta_outer)
  y_outer <- r2 * sin(theta_outer)
  
  # 2) Bottom right-half cap: connects (0, -r2) -> (0, -r1)
  phi_bot <- seq(-pi/2, pi/2, length.out = ncap)           # right half (x >= 0)
  x_bot_cap <- r_cap * cos(phi_bot)
  y_bot_cap <- -k + r_cap * sin(phi_bot)                # center at (0, -k)
  
  # 3) Inner left semicircle: bottom -> top (theta from 3pi/2 -> pi/2)
  theta_inner <- seq(3*pi/2, pi/2, length.out = nmain)
  x_inner <- r1 * cos(theta_inner)
  y_inner <- r1 * sin(theta_inner)
  
  # 4) Top right-half cap: connects (0, r1) -> (0, r2)
  phi_top <- seq(-pi/2, pi/2, length.out = ncap)           # right half (x >= 0)
  x_top_cap <- r_cap * cos(phi_top)
  y_top_cap <-  k + r_cap * sin(phi_top)                # center at (0, k)
  
  # Combine in order to form the closed loop
  x <- c(x_outer, x_bot_cap, x_inner, x_top_cap)
  y <- c(y_outer, y_bot_cap, y_inner, y_top_cap)
  
  t0 <- seq(0,1,length.out=length(x))
  cbind(x,y, cos(4*pi*t0)) + rnorm(3*length(x), sd=noise)
}

plot_datafit_1d3D_multi_static <- function(data, fitlist,theta=40,phi=20, periodic=F) {
  library(plot3D)
  
  S <- seq(0, 1, length.out = 1000)
  #cols <- rainbow(length(fitlist$spline_list))
  cols <- wacolors$rainier
  
  # Base scatter plot (data points)
  scatter3D(
    x = data[, 1], y = data[, 2],z = data[, 3],
    col = rgb(0, 0, 0, alpha=0.3), pch = 16, cex = 0.5,
    xlab = "X",
    ylab = "Y",
    zlab = "Z",
    theta = theta,
    phi = phi,
    ticktype = "detailed",
    box=TRUE,shade=0.8, bty="g"
  )
  
  # Add fitted spline curves
  last <- length(fitlist$spline_list)
  
  for (i in seq_along(fitlist$spline_list)) {
    fit <- fitlist$spline_list[[i]]
    splinevis <- fit(S)
    
    lines3D(
      x = splinevis[, 1],
      y = splinevis[, 2],
      z = splinevis[, 3],
      col = cols[i],
      lwd = ifelse((i ==last) & periodic, 7, 2),
      add = TRUE
    )
  }
  
  # Legend
  legend(
    "topright",
    legend = paste("λ =", fitlist$lambda_list),
    col = cols,
    lwd = 1,
    bty = "n"
  )
}


# Keep the historical seed and sampling sequence; discard the unrelated
# 2D surface experiments formerly preceding this section.
set.seed(1)
X1 <- sine1d3D(n_sample, rotation = pi/6, period = 2*pi, noise = 0.02)$X
result1 <- CompactPME::pme_1d_interval(
  X1, lambdas, parallel = FALSE, k_init = 5L, max_iter = fit_iterations,
  optimize_lambda = FALSE)
X2 <- helix1d3D(n_sample, rotation = pi/6, zscale = 3, noise = 0.02)$X
result2 <- CompactPME::pme_1d_interval(
  X2, lambdas, parallel = FALSE, k_init = 5L, max_iter = fit_iterations,
  optimize_lambda = FALSE)
png(file.path(figure_dir, "1d3D_sine_spiral.png"), width = 4800,
    height = 3000, res = 400)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 2))
plot_datafit_1d3D_multi_static(X1, result1, theta = 10, phi = 30)
plot_datafit_1d3D_multi_static(X2, result2, theta = 10, phi = 10)
dev.off()

X1 <- flower1d3D(n_sample, petals = 5, noise = 0.02)$X
result1 <- CompactPME::pme_1d_periodic(
  X1, lambdas, parallel = FALSE, max_iter = fit_iterations,
  optimize_lambda = FALSE)
X2 <- moon1d3D(n_sample, noise = 0.02)
result2 <- CompactPME::pme_1d_periodic(
  X2, lambdas, parallel = FALSE, max_iter = fit_iterations,
  optimize_lambda = FALSE)
png(file.path(figure_dir, "1d3D_star_moon.png"), width = 4800,
    height = 3000, res = 400)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 2))
plot_datafit_1d3D_multi_static(X1, result1, theta = 40, phi = 70, periodic = TRUE)
plot_datafit_1d3D_multi_static(X2, result2, theta = -70, phi = 40, periodic = TRUE)
dev.off()
