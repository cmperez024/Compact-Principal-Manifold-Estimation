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
library(splines)

#library(rgl)
#dos2unix
#squeue --start --me

library(CompactPME)
fibonaccisphere <- getFromNamespace("fibonaccisphere", "CompactPME")
spline_error <- getFromNamespace("spline_error", "CompactPME")
var_het <- getFromNamespace("var_het", "CompactPME")
local_var_1d <- getFromNamespace("local_var_1d", "CompactPME")
project_circle <- getFromNamespace("project_circle", "CompactPME")
project_grid <- getFromNamespace("project_grid", "CompactPME")
project_optimize <- getFromNamespace("project_optimize", "CompactPME")
iso_try <- getFromNamespace("iso_try", "CompactPME")
normalize <- getFromNamespace("normalize", "CompactPME")
Rk <- getFromNamespace("Rk", "CompactPME")
spline1d_interval <- getFromNamespace("spline1d_interval", "CompactPME")
spline1d_interval_approx <- getFromNamespace("spline1d_interval_approx", "CompactPME")
pmeClosed_1d_periodic <- CompactPME::pme_1d_periodic
krt <- getFromNamespace("krt", "CompactPME")
krst <- getFromNamespace("krst", "CompactPME")
spline1d_periodic <- getFromNamespace("spline1d_periodic", "CompactPME")
spline1d_periodic_approx <- getFromNamespace("spline1d_periodic_approx", "CompactPME")
pmeClosed_1d_interval <- CompactPME::pme_1d_interval

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
# =============== 1d Manifolds =============== 

# General functions (applies to both templates)

# Compute residuals

# Helper function for variance heterogeneity given smoothed residuals from local_var_1d


#  Compute smoothed residuals using kernel smoothing

# Circular projection, for initialization puprposes

# Optimization over a grid. This is to find the best projection index on a 1d Manifold

# Find projection index using an optimization routine. Grid method is preferred for global solutions

# An isomap method to find the smallest k so that the graph is not fragmented



# SPLINES AND PME: INTERVAL ============================

# The reproducing kernel for the interval

# The exact spline on the interval

# The approximate spline on the interval



# PME For the INTERVAL
# data should be NxD where D is 2 or 3.
# gridSize determines the partition for the projection grid method
# lvp_grid determines how many smoothed residual estimates for kernel regression; used for variance heterogeneity
# cl is a cluster, used for simulations
# parallel is for whether you want parallel computing or not (recommended for many lambda and larger N)

# PME for the CIRCLE
# ==================== Kernel functions ====================
# Bernoulli polynomial helper
# Bernoulli kernel (use this)


# Spline helper
# ========== Spline fitting ============================
# Given an orderd data set X in nxD (ordered with t) and penalty, fit the closed 1d spline
# Returns function f mapping [0,1] to R^D. Can input multiple t (as a vector) to evaluate each t value


# The approximate spline on the interval

# ========== Error for fit ====================

# PME method
# Parameters mostly identical to above.



#========== datasets ============
flower1d2D <- function(n, petals, noise, t){
  t0 <- if(missing(t)) sort(runif(n)) else t
  
  r <- (1 + 0.3 * sin(petals *2*pi* t0))
  return(list(t=t0, X = cbind(r * cos(2*pi*t0),  r * sin(2*pi*t0))+rnorm(2*n, sd=noise)))
}


#chatgpt
spiral1d2D_arclength <- function(n, period = 4*pi, scale = 1, noise = 0){
  
  omega <- period
  
  # arc length function
  s_fun <- function(t){
    0.5 * t * sqrt(1 + omega^2 * t^2) +
      asinh(omega * t) / (2 * omega)
  }
  
  # total length
  L <- s_fun(1)
  
  # sample uniformly in arc length
  s_star <- sort(runif(n, 0, L))
  
  # invert numerically
  t_star <- sapply(s_star, function(s_target){
    uniroot(function(t) s_fun(t) - s_target,
            interval = c(0,1))$root
  })
  
  # evaluate curve
  X <- cbind(scale * t_star * cos(omega * t_star),
             scale * t_star * sin(omega * t_star))
  
  # add noise
  if(noise > 0){
    X <- X + matrix(rnorm(2*n, sd=noise), ncol=2)
  }
  
  list(X = X, t = t_star)
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


halfflower1d2D <- function(n, petals, noise, t, rotation=0){
  t0 <- if(missing(t)) sort(runif(n)) else t
  
  r <- (1 + 0.3 * sin(petals *pi* t0))
  return(list(t=t0, X = cbind(r * cos(pi*t0),  r * sin(pi*t0)) %*% rotation2(rotation) +rnorm(2*n, sd=noise)))
} 

# ============================ RUN SIMULATION =================================
# Deliverables:
# A: TDA plot + fits + functional
# B: CV plot 
# 1d curve in 2D



# =============== 1d ======================
# 1. load data
set.seed(1)
N <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 30L else 2500L
boot <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 2L else 100L
lambda_N <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 3L else 150L
max_it <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 2L else 50L
min_it <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 1L else 10L


#dataset_interval <- spiral1d2D_arclength(N, period=5*pi, noise=0.02)
dataset_interval <- halfflower1d2D(N, petals=0, rotation=3*pi/2, noise=0.02)
dataset_periodic <- flower1d2D(N, petals=5, noise=0.02)

# 2. perform tda and save result
tda_interval <- runTDA_1d2D(dataset_interval$X, boot_size=boot)
tda_periodic <- runTDA_1d2D(dataset_periodic$X, boot_size=boot)
saveRDS(tda_interval, "results_1d/1d_interval_tda_local_exact.rds")
saveRDS(tda_periodic, "results_1d/1d_periodic_tda_local_exact.rds")

# 3. perform pme and save result
lambdas_int <- 10^(seq(-12, 3, length.out=lambda_N))
lambdas_per <- 10^(seq(-12, 3, length.out=lambda_N))
pme_interval <- pmeClosed_1d_interval(dataset_interval$X, lambdas_int, err_tol = 0.001, max_iter = max_it, min_iter = min_it, optimize_lambda=T,
                                      parallel = !identical(Sys.getenv("PME_SMOKE_TEST"), "1"), gridSize=N,approx=F, lvp_grid=10*N)

pme_periodic <- pmeClosed_1d_periodic(dataset_periodic$X, lambdas_per, err_tol = 0.001, max_iter = max_it, min_iter = min_it, optimize_lambda=T,
                                      parallel = !identical(Sys.getenv("PME_SMOKE_TEST"), "1"), gridSize=N,approx=F, lvp_grid=10*N)

saveRDS(pme_interval, "results_1d/1d_interval_pme_local_exact.rds")
saveRDS(pme_periodic, "results_1d/1d_periodic_pme_local_exact.rds")
