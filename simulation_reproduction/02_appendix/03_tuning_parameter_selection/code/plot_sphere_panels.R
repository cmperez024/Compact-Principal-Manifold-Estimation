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

# Regenerate this figure's own inputs when its cache is absent.
reproduce_inputs <- function() {
library(vegan) # isomap
library(plotly) # plotting
library(foreach) # parallelization
library(doParallel)
library(CompactPME)
krt <- CompactPME:::krt
krst <- CompactPME:::krst
Rk <- CompactPME:::Rk
qm <- CompactPME:::qm
expr <- CompactPME:::expr
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
N <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 30L else 2500L
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
lambda_N <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 3L else 100L
lambdas <- 10^(seq(-14,-3,length.out=lambda_N))
pme_sphere <- pme_2d_periodic(dataset$X, lambdas, optimize_lambda = T,parallel = !identical(Sys.getenv("PME_SMOKE_TEST"), "1"), max_iter = if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 2L else 50L,
  err_tol = if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 0.01 else 1e-8,
  resid_smooth_size = 5*N)
pmepath <- paste0("results_sph/flower_pme", job_label, ".rds")
saveRDS(pme_sphere, pmepath)


N <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 30L else 2500L
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
lambda_N <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 3L else 100L
lambdas <- 10^(seq(-14,-3,length.out=lambda_N))
pme_sphere <- pme_2d_periodic(dataset$X, lambdas, optimize_lambda = T,parallel = !identical(Sys.getenv("PME_SMOKE_TEST"), "1"),max_iter = if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 2L else 50L,
  err_tol = if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 0.01 else 1e-8,
  resid_smooth_size = 5*N)
pmepath <- paste0("results_sph/cashew_pme", job_label, ".rds")
saveRDS(pme_sphere, pmepath)


}
needed_inputs <- file.path("results_sph", c("flower_tda.rds", "flower_pme.rds", "cashew_tda.rds", "cashew_pme.rds"))
if (!all(file.exists(needed_inputs))) reproduce_inputs()

library(pracma)
library(dplyr)
library(ggplot2)
library(TDA)
library(plot3D)
library(tidyr)
library(patchwork)
library(png)
library(grid)
library(gridExtra)

# 2d TDA
set.seed(1)
fibonaccisphere = function(n=1000, r=1, out.xyz=TRUE, out.sph=FALSE) {

  if (n<1 | round(n)!=n) stop('n must be a positive integer')
  if (!out.xyz & !out.sph) stop('either out.xyz and/or out.sph must be TRUE')

  goldenratio = (1+sqrt(5))/2
  i = seq(n)-0.5
  z = 1-2*i/n # z-coordinate for unit sphere
  theta = acos(pmax(-1,pmin(1,z))) # polar angle
  phi = (2*pi*i/goldenratio)%%(2*pi)

  if (out.xyz) {
    x = r*sin(theta)*cos(phi)
    y = r*sin(theta)*sin(phi)
    z = r*z
    if (out.sph) {
      out = cbind(x=x,y=y,z=z,theta=theta,phi=phi)
    } else {
      out = cbind(x=x,y=y,z=z)
    }
  } else {
    out = cbind(theta=theta,phi=phi)
  }

  return(out)

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

tda_2_plot2d <- function(tda_res, label="", hide_legend=F, title=""){
  Diag <- tda_res$Diag
  boot <- tda_res$boot
  alpha <- tda_res$alpha

  conf_label <- paste0(round((1 - alpha) * 100), "% confidence band")

  pd_df <- Diag$diagram %>% unclass() %>% as.data.frame()

  grid_points <- seq(min(pd_df$Birth)-1, max(pd_df$Birth)+5, length.out = 100)

  #ribbon_df <- data.frame(
  #  x = c(grid_points - boot, rev(grid_points + boot)),
  #  y = c(grid_points + boot, rev(grid_points - boot)),
  #  Band=factor(conf_label))

  ribbon_df <- data.frame(
    x= c(grid_points, rev(grid_points)),
    y=c(grid_points + 2*boot,rev(grid_points)),
    Band = conf_label
  )


  xmax <- max(pd_df$Birth)*1.05
  ymax <- max(pd_df$Death)*1.05

  legendpos <- ifelse(hide_legend, "none", "bottom")

  pd_df %>%
    ggplot() +
    geom_point(aes(x=Birth,y=Death, color = factor(dimension), shape=factor(dimension)),size=2) +
    scale_color_manual(name=paste0(label, "Dimension"),values=c("0"="black", "1"="red", "2"= "blue")) +
    scale_shape_manual(name=paste0(label, "Dimension"),values = c("0" = 16, "1" = 2, "2"= 5)) +
    geom_abline(slope=1,intercept=0)  +
    geom_polygon(data=ribbon_df, aes(x=x,y=y, fill = Band),alpha=0.2) +
    theme_bw() +
    theme(panel.grid=element_blank(), legend.position=legendpos,legend.background=element_rect(fill="white",color="black"))+
    scale_fill_manual(values = setNames("red", levels(ribbon_df$Band)), name = paste0(label, "Band"))+theme(aspect.ratio=1)+
    #xlim(c(0,xmax))+ylim(c(0,ymax))
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) + ggtitle(title)
}


flower2d3D_unif <- function(N = 300,
                            r0 = 1, a = 0.3, petals = 5,
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


flower2d3D_func <- function(t,
                            r0 = 1, a = 0.3, petals = 6,
                            b = 0.5
) {

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

  cbind(x, y, z)
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

unifmoon2d3D_func <- function(t, bend = pi,
                         radius = 1, bend_radius = 1.2, noise=0) {

  S <- t/ sqrt(rowSums(t^2))

  x <- radius * S[,1]
  y <- radius * S[,2]
  z <- radius * S[,3]

  # 2. bend map
  new_x <- cos(bend * x / (2*radius)) * (bend_radius * radius + y)
  new_y <- sin(bend * x / (2*radius)) * (bend_radius * radius + y)
  new_z <- z

  cbind(new_x, new_y, new_z)
}
# 2d Viz plots

genSurface <- function(f){
  M <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 10L else 50L
  thetas <- seq(from = 0, to = 2*pi, length.out=M)
  phis <- seq(from =-pi/2, to = pi/2, length.out=M)
  rr <- 1
  grid <- expand.grid(thetas,phis,1)
  cartgrid <- t(apply(grid, 1, sph2cart))

  # Evaluate manifold
  fvalues <- f(cartgrid)

  Xmat <- matrix(fvalues[,1], nrow = M, ncol = M, byrow = TRUE)
  Ymat <- matrix(fvalues[,2], nrow = M, ncol = M, byrow = TRUE)
  Zmat <- matrix(fvalues[,3], nrow = M, ncol = M, byrow = TRUE)

  list(X=Xmat,Y=Ymat,Z=Zmat)
}

plot_datafit_2d3D_static <- function(X, fit, theta=40,phi=20, bty="n", title="Data and Surface",
                                     point_size = 0.5){
  surf <- genSurface(fit)



  # Get the lim given by min over all dimensions and max over all dimensions
  #lims <- c(min(extract_lims[1,]), max(extract_lims[2,]))
  lims <- c(min(X), max(X))

  # Start with data points
  scatter3D(
    X[,1], X[,2], X[,3],
    pch = 16, cex = point_size, col = "black",
    theta = theta, phi = phi,
    main = title,
    bty = bty,
    ticktype = "detailed", xlim=lims, ylim=lims,zlim=lims


  )

  # Overlay surface
  surf3D(
    surf$X, surf$Y, surf$Z,
    add = TRUE,
    colvar = surf$Z, # overlay on existing plot
    # col = "Viridis",
    alpha = 0.6,             # semi-transparent
    border = "grey",
    #xlim=c(-5,5),
    bty = bty,
    colkey=F,

    xlim=lims, ylim=lims,zlim=lims

  )


}

# choose lambda
# get plots of three choices of lambda
get_closest_lambda_ind <- function(result, lambda_target) {

  idx <- sapply(lambda_target, function(l) {
    which.min(abs(result$lambda_list - l))
  })

  unique(idx)
}


## CV Plots


# Plot all kernels or just pelletier?
cv_plots2d <- function(result, f_true, cutoff= 5e-4){


  reshapedf <- function(df){
    as.data.frame(df)  %>% pivot_longer(cols=-lambda,names_to="kernel_method",values_to="value")
  }

    # Results for pelletier only

    meanplt <- result$mean %>% reshapedf() %>% filter(kernel_method=="pelletier") %>% ggplot(aes(x=lambda, y=value)) + geom_point()+ geom_line()+ scale_x_log10() + theme_bw()+
      theme(panel.grid=element_blank()) + ylab("Volume-averaged conditional squared residual") + xlab("λ")+ ggtitle("(ii)") + geom_vline(xintercept=cutoff, linetype="dashed")


    sdplt <- result$sd %>% reshapedf() %>% filter(kernel_method=="pelletier") %>% ggplot(aes(x=lambda, y=value)) + geom_point()+ geom_line()+ scale_x_log10() + theme_bw()+
      theme(panel.grid=element_blank()) + ylab("Standard Deviation") + xlab("λ")+ ggtitle("(iii)")

    # Min_df
    mindf <- result$cv %>% reshapedf() %>% filter(kernel_method=="pelletier") %>% group_by(kernel_method) %>% slice_min(value, with_ties=F)
    # thresholded min
    mindf2 <- result$cv %>% reshapedf() %>% filter(lambda < cutoff, kernel_method=="pelletier") %>% group_by(kernel_method) %>% slice_min(value, with_ties=F)

    mindfs <- rbind(cbind(mindf, optimal="naive"), cbind(mindf2, optimal="thresholded"))


    cvplt <-result$cv %>% reshapedf() %>%  filter(kernel_method=="pelletier") %>%
      ggplot(aes(x=lambda, y=value)) + geom_point()+ geom_line()+ scale_x_log10() +
      theme_bw()+ geom_vline(xintercept = mindf2$lambda, col="red",linetype="dashed")+
      theme(panel.grid=element_blank(), legend.position="none")+
      ylab("Coefficient of Variation") + xlab("λ") + ggtitle("(i)")



  #(cvplt + meanplt)/ (sdplt+ plot_spacer())
  list(cv = cvplt, mean = meanplt, sd = sdplt)
}


# Plots of surfaces with optimal lambdas ======
get_closest_lambda_ind <- function(result, lambda_target) {

  idx <- sapply(lambda_target, function(l) {
    which.min(abs(result$lambda_list - l))
  })

  idx
}

reshapedf <- function(df){
  as.data.frame(df)  %>% pivot_longer(cols=-lambda,names_to="kernel_method",values_to="value")
}





#=========== Start FLOWER =============
# Load / setup data
tda_sphere <- readRDS("results_sph/flower_tda.rds")
(p_tda_sph <- tda_2_plot2d(tda_sphere,title="(iv)",hide_legend=T))
#set.seed(1)
ggsave("panel_flower/panel_A_4.png", width=3200, height=3200, units = "px",dpi=600)


pme_sphere <- readRDS("results_sph/flower_pme.rds")

cutoff <- 1e-4
ftrue <- flower2d3D_unif(nrow(pme_sphere$dataset))
cv_plots2d(pme_sphere, ftrue$X, cutoff=cutoff)
#ggsave("panel_flower/cv_res.png", width=4000, height=3200, units = "px",dpi=300)
cvplots <- cv_plots2d(pme_sphere, ftrue$X, cutoff=cutoff)

cvplots$cv
ggsave("panel_flower/panel_A_1.png", width=3200, height=3200, units = "px",dpi=600)
cvplots$mean
ggsave("panel_flower/panel_A_2.png", width=3200, height=3200, units = "px",dpi=600)
cvplots$sd
ggsave("panel_flower/panel_A_3.png", width=3200, height=3200, units = "px",dpi=600)



# Min_df
mindf <- pme_sphere$cv %>% reshapedf() %>% filter(kernel_method=="pelletier") %>% group_by(kernel_method) %>% slice_min(value, with_ties=F)
# thresholded min
mindf2 <- pme_sphere$cv %>% reshapedf() %>% filter(lambda < cutoff, kernel_method=="pelletier") %>% group_by(kernel_method) %>% slice_min(value, with_ties=F)

mindfs <- rbind(cbind(mindf, type="naive"), cbind(mindf2, type="cutoff"))

lam_min <- mindfs$lambda
kern_min <- mindfs$type

lam_min_ind <- get_closest_lambda_ind(pme_sphere, lam_min)

ilabs <- c("(ii)", "(iii)")
for(i in 1:length(lam_min)){

  png(paste0("panel_flower/panel_B_",i,".png"), width = 3200, height=3200, res=600)
  par(mar=c(0.02,0.1,3,0.1))

  lam_i <- lam_min_ind[i]
  mytitle <- paste0(ilabs[i]," ", kern_min[i], ", λ = ",formatC(pme_sphere$lambda_list[lam_i],
                                                                format = "e", digits=2))
  plot_datafit_2d3D_static(pme_sphere$dataset, pme_sphere$spline_list[[lam_i]], phi=40, point_size=0.4,
                           title=mytitle)

  dev.off()
}

lams <- c(10^-14, mindf2$lambda ,10^-3)
lams_ind <- get_closest_lambda_ind(pme_sphere, lams)

plotz <- vector("list", length=length(lams))

prelab <- c("(vi)", "(vii)", "(viii)")

for(i in 1:length(lams)){

  png(paste0("panel_flower/panel_B_",i+1,".png"), width = 3200, height=3200, res=600)
  par(mar=c(0.02,0.1,3,0.1))
  lam_i <- lams_ind[i]
  plot_datafit_2d3D_static(pme_sphere$dataset, pme_sphere$spline_list[[lam_i]], phi=40, point_size=0.4,
                           title=paste0(prelab[i], " λ = ",formatC(pme_sphere$lambda_list[lam_i],
                                                            format = "e", digits=2)))

  dev.off()
}
png(paste0("panel_flower/panel_B_1.png"), width = 3200, height=3200, res=600)
ftrue <- function(t) flower2d3D_func(t, petals=5)
par(mar=c(0.02,0.1,3,0.1))
plot_datafit_2d3D_static(pme_sphere$dataset, ftrue, phi=40, point_size=0.4, title="(v) Data and True Surface")
dev.off()
#plot_datafit_2d3D(pme_sphere$dataset, ftrue)



filesA <- sprintf("panel_flower/panel_A_%d.png", 1:4)
filesB <- sprintf("panel_flower/panel_B_%d.png", 1:4)
files <- c(filesA, filesB)
# Read all images
imgs <- lapply(files, function(f) {
  rasterGrob(readPNG(f), interpolate = TRUE)
})
# Arrange in a 2x2 grid
png(file.path(figure_dir, "flower_panel.png"), width = 12800, height = 6400,res=600)
grid.arrange(grobs = imgs, ncol = 4, nrow = 2)
dev.off()


# === Start CASHEW ==============

# Load / setup data
tda_sphere <- readRDS("results_sph/cashew_tda.rds")
(p_tda_sph <- tda_2_plot2d(tda_sphere,title="(iv)",hide_legend=T))
#set.seed(1)
ggsave("panel_cashew/panel_A_4.png", width=3200, height=3200, units = "px",dpi=600)


pme_sphere <- readRDS("results_sph/cashew_pme.rds")

cutoff <- 5e-4
#ftrue <- flower2d3D_unif(1000)
ftrue <- unifmoon2d3D(N=nrow(pme_sphere$dataset))
#cv_plots2d(pme_sphere, ftrue$X, cutoff=cutoff)
#ggsave("panel_cashew/cv_res.png", width=4000, height=3200, units = "px",dpi=300)
cvplots <- cv_plots2d(pme_sphere, ftrue$X, cutoff=cutoff)

cvplots$cv
ggsave("panel_cashew/panel_A_1.png", width=3200, height=3200, units = "px",dpi=600)
cvplots$mean
ggsave("panel_cashew/panel_A_2.png", width=3200, height=3200, units = "px",dpi=600)
cvplots$sd
ggsave("panel_cashew/panel_A_3.png", width=3200, height=3200, units = "px",dpi=600)



# Min_df
mindf <- pme_sphere$cv %>% reshapedf() %>% filter(kernel_method=="pelletier") %>% group_by(kernel_method) %>% slice_min(value, with_ties=F)
# thresholded min
mindf2 <- pme_sphere$cv %>% reshapedf() %>% filter(lambda < cutoff, kernel_method=="pelletier") %>% group_by(kernel_method) %>% slice_min(value, with_ties=F)

mindfs <- rbind(cbind(mindf, type="naive"), cbind(mindf2, type="cutoff"))

lam_min <- mindfs$lambda
kern_min <- mindfs$type

lam_min_ind <- get_closest_lambda_ind(pme_sphere, lam_min)

ilabs <- c("(ii)", "(iii)")
for(i in 1:length(lam_min)){

  png(paste0("panel_cashew/panel_B_",i,".png"), width = 3200, height=3200, res=600)
  par(mar=c(0.02,0.1,3,0.1))

  lam_i <- lam_min_ind[i]
  mytitle <- paste0(ilabs[i]," ", kern_min[i], ", λ = ",formatC(pme_sphere$lambda_list[lam_i],
                                                                format = "e", digits=2))
  plot_datafit_2d3D_static(pme_sphere$dataset, pme_sphere$spline_list[[lam_i]], phi=40, point_size=0.4,
                           title=mytitle)

  dev.off()
}

lams <- c(10^-14, mindf2$lambda ,10^-3)
lams_ind <- get_closest_lambda_ind(pme_sphere, lams)

plotz <- vector("list", length=length(lams))

prelab <- c("(vi)", "(vii)", "(viii)")

for(i in 1:length(lams)){

  png(paste0("panel_cashew/panel_B_",i+1,".png"), width = 3200, height=3200, res=600)
  par(mar=c(0.02,0.1,3,0.1))
  lam_i <- lams_ind[i]
  plot_datafit_2d3D_static(pme_sphere$dataset, pme_sphere$spline_list[[lam_i]], phi=40, point_size=0.4,
                           title=paste0(prelab[i], " λ = ",formatC(pme_sphere$lambda_list[lam_i],
                                                                   format = "e", digits=2)))

  dev.off()
}
png(paste0("panel_cashew/panel_B_1.png"), width = 3200, height=3200, res=600)
ftrue <- function(t) unifmoon2d3D_func(t)
par(mar=c(0.02,0.1,3,0.1))
plot_datafit_2d3D_static(pme_sphere$dataset, ftrue, phi=40, point_size=0.4, title="(v) Data and True Surface")
dev.off()
#plot_datafit_2d3D(pme_sphere$dataset, ftrue)



filesA <- sprintf("panel_cashew/panel_A_%d.png", 1:4)
filesB <- sprintf("panel_cashew/panel_B_%d.png", 1:4)
files <- c(filesA, filesB)
# Read all images
imgs <- lapply(files, function(f) {
  rasterGrob(readPNG(f), interpolate = TRUE)
})
# Arrange in a 2x2 grid
png(file.path(figure_dir, "cashew_panel.png"), width = 12800, height = 6400,res=600)
grid.arrange(grobs = imgs, ncol = 4, nrow = 2)
dev.off()




