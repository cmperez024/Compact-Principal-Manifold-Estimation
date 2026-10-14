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

library(vegan) # isomap
library(plotly) # plotting
library(foreach) # parallelization
library(doParallel)
library(pracma)
library(TDA)
#library(rgl)
#dos2unix
#squeue --start --me


library(CompactPME)
fibonaccisphere <- getFromNamespace("fibonaccisphere", "CompactPME")
iso_try <- getFromNamespace("iso_try", "CompactPME")
var_het <- getFromNamespace("var_het", "CompactPME")
d_circ_vec <- getFromNamespace("d_circ_vec", "CompactPME")
pelletier_kernel <- getFromNamespace("pelletier_kernel", "CompactPME")
qm <- getFromNamespace("qm", "CompactPME")
expr <- getFromNamespace("expr", "CompactPME")
Rpp <- getFromNamespace("Rpp", "CompactPME")
spline_error <- getFromNamespace("spline_error", "CompactPME")
spline2d <- getFromNamespace("spline2d", "CompactPME")
project_optimize2 <- getFromNamespace("project_optimize2", "CompactPME")
pme_2d_periodic <- getFromNamespace("pme_2d_periodic", "CompactPME")

#==== HELPER FUNCTIONS ==================
#cooltools functions necessary. Taken from github to avoid loading unnecessary functions




rotation2 = function(angle) {
  c = cos(angle)
  s = sin(angle)
  return(rbind(c(c,-s),c(s,c)))
  
}

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

# =============== 2d Manifolds =============== 

# Some variance heterogeneity methods===============


# Local variance helper functions
local_var_2d <- function(Z, resid_sq, eval_size = 1000,bw_grid = 100) {
  xriem <- wrap.sphere(Z) # turn our dependent variable to riemann sphere object
  
  # bws <- seq(from = 0, to = 1, length.out= grid_eval)[-1]# get set of bandwidths to test
  bws <- seq(from = 0.1, to = 1, length.out= bw_grid) # bw too small results in NA
  
  intrin <- riem.m2skregCV(xriem, resid_sq, geometry="intrinsic", bandwidths =bws)
  extrin <- riem.m2skregCV(xriem, resid_sq, geometry="extrinsic",bandwidths =bws)
  
  # Predictions on new set
  xnew <- wrap.sphere(fibonaccisphere(eval_size)) # evenly spaced sphere points
  
  # Get smoothed results
  smooth_i <- predict(intrin, xnew)
  smooth_e <- predict(extrin, xnew)
  
  # return CV and bandwidth
  list(s2_intrin =  smooth_i, bw_intrin = intrin$bandwidth,
       s2_extrin = smooth_e, bw_extrin = extrin$bandwidth)
}


# Methods for pelletier kernel

# Gaussian kernel

#' Pelletier, B. (2006). Non-parametric regression estimation on closed Riemannian manifolds. Journal of Nonparametric Statistics, 18(1), 57-67.


# Splines and intialization functions==================
# Kernel functions
m <- 2

# qm function
# kernelhelper


# Kernel function

# error helper

# Spline fit helper

# Projection via optimization
# Optimize surface based on spline values and previous projection parameters
# (cartesian)

# PME for the 2-sphere =================

#========== datasets ============
flower1d2D <- function(n, petals, noise, t){
  t0 <- if(missing(t)) sort(runif(n)) else t
  
  r <- (1 + 0.3 * sin(petals *2*pi* t0))
  return(list(t=t0, X = cbind(r * cos(2*pi*t0),  r * sin(2*pi*t0))+rnorm(2*n, sd=noise)))
}

halfflower1d2D <- function(n, petals, noise, t){
  t0 <- if(missing(t)) sort(runif(n, max=0.5)) else t
  
  r <- (1 + 0.3 * sin(petals *2*pi* t0))
  return(list(t=t0, X = cbind(r * cos(2*pi*t0),  r * sin(2*pi*t0))+rnorm(2*n, sd=noise)))
} 

flower2d3D_unif <- function(N = 300,
                            r0 = 1, a = 0.3, petals = 6,
                            b = 0.5, noise = 0.01,
                            t = NULL) {
  
  # If no sphere provided, generate uniform Fibonacci sphere
  if (is.null(t)) {
    t <- fibonaccisphere(N)  # N x 3 matrix, uniform on S^2
  }
  
  # Ensure unit sphere
  t <- t / sqrt(rowSums(t^2))
  
  # Convert to spherical coordinates
  coords <- cart2sph(t)
  phi   <- coords[, "theta"]     # longitude
  elev  <- coords[, "phi"]       # elevation
  z     <- t[, 3]                # vertical coordinate
  
  # Radial modulation in xy-plane
  r_xy <- r0 * (1 + a * cos(petals * phi))
  
  # Flower deformation
  x <- r_xy * sqrt(1 - z^2) * cos(phi)
  y <- r_xy * sqrt(1 - z^2) * sin(phi)
  z <- b * z
  
  pts <- cbind(x, y, z) + matrix(rnorm(3*N, sd = noise), N, 3)
  
  return(list(t = t, X = pts))
}


unifmoon2d3D <- function(N, oversample = 5, bend = pi,
                         radius = 1, bend_radius = 1.2, noise=0) {
  
  # 1. oversample sphere uniformly
  S <- fibonaccisphere(N * oversample)
  
  x <- radius * S[,1]
  y <- radius * S[,2]
  z <- radius * S[,3]
  
  # 2. bend map
  new_x <- cos(bend * x / (2*radius)) * (bend_radius * radius + y)
  new_y <- sin(bend * x / (2*radius)) * (bend_radius * radius + y)
  new_z <- z
  
  pts <- cbind(new_x, new_y, new_z)
  
  # 3. Rejection thinning for uniform coverage
  # compute kernel density to approximate sampling density
  library(FNN)
  k <- 10
  d <- get.knn(pts, k)$nn.dist[,k]       # k-th NN distance
  w <- d^3                               # local volume ~ r^3
  prob <- w / max(w)
  
  keep <- sample(seq_len(nrow(pts)), N, prob=prob)
  
  list(X = pts[keep,] + rnorm(3*N, sd=noise),
       t = S[keep,])
}


#======tda========
runTDA_1d2D <- function(data, alpha=0.01, boot_size=30, by = 0.01)
{
  Xlim <- c(min(data[,1]), max(data[,1]))
  Ylim <- c(min(data[,2]), max(data[,2]))
  lim <- cbind(Xlim, Ylim)
  
  Diag <- gridDiag(data, distFct, lim = lim, by = by)
  boot <- bootstrapDiagram(X=data, FUN=distFct, lim = lim,by=by, parallel = !identical(Sys.getenv("PME_SMOKE_TEST"), "1"),alpha=alpha,B=boot_size)
  
  list(Diag=Diag, boot=boot, alpha =alpha)
}

runTDA_2d3D <- function(data, alpha=0.01, boot_size=30)
{
  Xlim <- c(min(data[,1]), max(data[,1]))
  Ylim <- c(min(data[,2]), max(data[,2]))
  Zlim <- c(min(data[,3]), max(data[,3]))
  lim <- cbind(Xlim, Ylim,Zlim)
  by <- 0.05
  
  Diag <- gridDiag(data, distFct, lim = lim, by = by)
  boot <- bootstrapDiagram(X=data, FUN=distFct, lim = lim,by=by, parallel = !identical(Sys.getenv("PME_SMOKE_TEST"), "1"),alpha=alpha,B=boot_size)
  
  list(Diag=Diag, boot=boot, alpha =alpha)
}




# ============================ RUN SIMULATION =================================
# Deliverables:
# A: TDA plot + fits + functional
# B: CV plot 
# 1d curve in 2D



# ======= 2d ==========
#rm(list = ls())

job_label <- ""

# 1. Load dataset
set.seed(1)
N <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 30L else 250L
dataset <- flower2d3D_unif(N, petals=5, noise=0.02)
#plot_data_3D(dataset$X)

# 2. TDA step
# Restore flower persistence calculation needed by the panel.
boot <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 2L else 100L
tda_sphere <- runTDA_2d3D(dataset$X, boot_size=boot)
#saveRDS()
#tda_2_plot_3D(tdares)
tdapath <- paste0("results_sph/flower_tda", job_label, ".rds")
saveRDS(tda_sphere, tdapath)

# 3. Run PME
lambda_N <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 3L else 10L
lambdas <- 10^(seq(-14,-3,length.out=lambda_N))
pme_sphere <- pme_2d_periodic(dataset$X, lambdas, optimize_lambda = T,parallel = !identical(Sys.getenv("PME_SMOKE_TEST"), "1"), resid_smooth_size = 5*N)
pmepath <- paste0("results_sph/flower_pme", job_label, ".rds")
saveRDS(pme_sphere, pmepath)


N <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 30L else 250L
dataset <- unifmoon2d3D(N,noise=0.02)
#plot_data_3D(dataset$X)

# 2. TDA step
boot <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 2L else 10L
tda_sphere <- runTDA_2d3D(dataset$X, boot_size=boot)
#saveRDS()
#tda_2_plot_3D(tdares)
tdapath <- paste0("results_sph/cashew_tda", job_label, ".rds")
saveRDS(tda_sphere, tdapath)

# 3. Run PME
lambda_N <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 3L else 10L
lambdas <- 10^(seq(-14,-3,length.out=lambda_N))
pme_sphere <- pme_2d_periodic(dataset$X, lambdas, optimize_lambda = T,parallel = !identical(Sys.getenv("PME_SMOKE_TEST"), "1"),resid_smooth_size = 5*N)
pmepath <- paste0("results_sph/cashew_pme", job_label, ".rds")
saveRDS(pme_sphere, pmepath)

