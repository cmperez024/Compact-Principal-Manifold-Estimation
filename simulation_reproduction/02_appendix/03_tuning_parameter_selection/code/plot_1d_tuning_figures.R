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
}
needed_inputs <- file.path("results_1d", c("1d_interval_tda_local_exact.rds", "1d_periodic_tda_local_exact.rds", "1d_interval_pme_local_exact.rds", "1d_periodic_pme_local_exact.rds"))
if (!all(file.exists(needed_inputs))) reproduce_inputs()

library(pracma)
library(dplyr)
library(purrr)
library(ggplot2)
library(TDA)
library(RColorBrewer)
library(patchwork)




#========== datasets ============
flower1d2D <- function(n, petals, noise, t){
  t0 <- if(missing(t)) sort(runif(n)) else t

  r <- (1 + 0.3 * sin(petals *2*pi* t0))
  return(list(t=t0, X = cbind(r * cos(2*pi*t0),  r * sin(2*pi*t0))+rnorm(2*n, sd=noise)))
}

halfflower1d2D <- function(n, petals, noise, t, rotation=0){
  t0 <- if(missing(t)) sort(runif(n)) else t

  r <- (1 + 0.3 * sin(petals *pi* t0))
  return(list(t=t0, X = cbind(r * cos(pi*t0),  r * sin(pi*t0)) %*% rotation2(rotation) +rnorm(2*n, sd=noise)))
}


rotation2 = function(angle) {
  c = cos(angle)
  s = sin(angle)
  return(rbind(c(c,-s),c(s,c)))

}

# Section 3 Plots

## Figure 3.2

### TDA

runTDA_1d2D <- function(data, alpha=0.01, boot_size=30, use_kde=F, by = 0.01)
{
  Xlim <- c(min(data[,1]), max(data[,1]))
  Ylim <- c(min(data[,2]), max(data[,2]))
  lim <- cbind(Xlim, Ylim)


  Diag <- gridDiag(data, distFct, lim = lim, by = by)
  boot <- bootstrapDiagram(X=data, FUN=distFct, lim = lim,by=by, parallel = !identical(Sys.getenv("PME_SMOKE_TEST"), "1"),alpha=alpha,B=boot_size)


  list(Diag=Diag, boot=boot, alpha =alpha)
}

tda_2_plot <- function(tda_res, label="", hide_legend=F, title=""){
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
    scale_color_manual(name=paste0(label, "Dimension"),values=c("0"="black", "1"="red")) +
    scale_shape_manual(name=paste0(label, "Dimension"),values = c("0" = 16, "1" = 2)) +
    geom_abline(slope=1,intercept=0)  +
    geom_polygon(data=ribbon_df, aes(x=x,y=y, fill = Band),alpha=0.2) +
    theme_bw() +
    theme(panel.grid=element_blank(), legend.position=legendpos,legend.background=element_rect(fill="white",color="black"))+
    scale_fill_manual(values = setNames("red", levels(ribbon_df$Band)), name = paste0(label, "Band"))+theme(aspect.ratio=1)+
    #xlim(c(0,xmax))+ylim(c(0,ymax))
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) + ggtitle(title)
}






tda_interval <- readRDS("results_1d/1d_interval_tda_local_exact.rds")
tda_periodic <- readRDS("results_1d/1d_periodic_tda_local_exact.rds")


# Get the plots
(p_tda_int <- tda_2_plot(tda_interval, title="(a)", label="(a, d) "))
#plot(tda_interval$Diag$diagram, band=2*tda_interval$boot)
(p_tda_per <- tda_2_plot(tda_periodic, title="(d)", label="(a, d) "))





### Manifold fit

pme_to_df <- function(result, param){
  # Return iteration details=======
  df_list <- Map(function(lambda, values) {

    data.frame(
      lambda    = lambda,
      iteration = seq_len(length(values)),
      cost_functional     = values
    )
  },
  result$lambda_list,
  result$cost_functional_list)

  iteration <- do.call(rbind, df_list)

  # return spline fit details=======
  S <- NULL
  if(param=="interval"){
    S <- seq(0,1,length.out=1000)
  }else if(param=="S1"){
    S <- seq(0,0.999,length.out=1000)
  }else{
    S <- fibonaccisphere(1000)
  }
  D <- ncol(result$spline_list[[1]](S))
  df <- do.call(rbind, lapply(seq_along(result$spline_list), function(i) {

    f      <- result$spline_list[[i]]
    lambda <- result$lambda_list[i]

    vals <- f(S)   # n × d matrix
    colnames(vals) <- c("x", "y", if (D == 3) "z")[seq_len(D)]

    data.frame(
      lambda = lambda,
      t      = S,
      vals
    )
  }))
  list(iteration = iteration, image = df)
}

# provide a list of target lambda and provide closest match (useful when want to plot
# select lambda in a fine grid of lambda)
get_closest_lambda_ind <- function(result, lambda_target) {

  idx <- sapply(lambda_target, function(l) {
    which.min(abs(result$lambda_list - l))
  })

  unique(idx)
}

# last_point_threshold draws the fitted lines with corresponding lambda above threshold to visualize collapse to mean
image_plot <- function(result, param, lambda_factor=T,lambda_target=NULL, last_point_threshold =Inf, sci_dec = 2,
                       title="", label=""){
  if(is.null(lambda_target)){
    lambda_target <- result$lambda_list
  }
  df <- pme_to_df(result, param)$image
  fact <- ifelse(lambda_factor, as.factor, function(x){x})


  lambda_subset <-result$lambda_list[get_closest_lambda_ind(result, lambda_target)]

  filtered_df <- df %>% filter(lambda %in% lambda_subset)
  point_df <- df %>% filter(lambda %in% lambda_subset & lambda >= last_point_threshold)



  lab <- paste0(label,"λ")
  data.frame(result$dataset) %>% `colnames<-`(c("x","y")) %>% ggplot(aes(x=x,y=y)) +
    geom_point(alpha=0.15) + geom_path(data=filtered_df, aes(x=x,y=y,col=fact(lambda)), size=0.85) + theme_bw() +
    xlab("x") + ylab("y")  +
    geom_point(data=point_df,aes(x=x,y=y, col=as.factor(lambda)),show.legend=F) +
    scale_color_brewer(palette ="Dark2", name =lab, labels = function(x) formatC(as.numeric(x),
                                                                                 format = "e",
                                                                                 digits = sci_dec))+
    theme(legend.position=c(0.98,0.98), legend.justification = c("right","top"),
          legend.background=element_rect(fill="white",color="black"), panel.grid=element_blank(),aspect.ratio=1)+
    ggtitle(title)
}




pme_interval <- readRDS("results_1d/1d_interval_pme_local_exact.rds")
pme_periodic <- readRDS("results_1d/1d_periodic_pme_local_exact.rds")



lams <- 10^c(-12, -7, 2)
mylab <- "(b, e) "
(p_pme_int <- pme_interval %>% image_plot(param="interval",lambda_target =lams, title="(b)", label=mylab))
(p_pme_per <-pme_periodic %>% image_plot(param="S1",lambda_target =lams, last_point_threshold = 1e-1,
                                         title="(e)",label=mylab))




### Cost function


result <- pme_periodic
param = "periodic"

library(scales) # for pretty breaks
iteration_plot <- function(result, param, lambda_optimal=NULL, title=""){
  if(is.null(lambda_optimal)){
    # lambda_optimal <- result$optimal
  }


  df <- pme_to_df(result, param)$iteration
  #fact <- ifelse(lambda_factor, as.factor, function(x){x})

  # Grab the optimal lambda
  lambda_subset <- result$lambda_list[get_closest_lambda_ind(result, lambda_optimal)]
  #print(lambda_subset)

  df_new <- df %>% filter(lambda==lambda_subset)

  df_new %>% ggplot(aes(x=iteration, y = cost_functional)) + geom_point() + geom_line() + theme_bw()+ theme(panel.grid=element_blank())+ ggtitle(title) +
    #scale_x_continuous(breaks=breaks_pretty(max(df_new %>% pull(iteration) %>% max()))) +
    xlab("Iteration") + ylab("Cost Functional")+ ggtitle(paste0(title, " λ = ", formatC(lambda_subset, format="e", digits=2)) )+
    theme(panel.grid=element_blank(),aspect.ratio=1) #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# plots===================
(p_func_int <- pme_interval %>% iteration_plot(param="interval",lambda_optimal= 1e-7, title="(c)"))

(p_func_per <- pme_periodic %>% iteration_plot(param="S1",lambda_optimal= 1e-7, title="(f)"))


### Combine and save

(plt<-(p_tda_int + p_pme_int + p_func_int)/(p_tda_per + p_pme_per + p_func_per)+ plot_layout(guides = 'collect') & theme(legend.position="bottom"))
ggsave(file.path(figure_dir, "interval_S1_full.png"), plt, width=9, height=6, dpi=800)





## Figure 3.3




# Section 5 Plots

## Figure 5.1

#stackoverflow
# Source - https://stackoverflow.com/a/34205938
# Posted by stas g
# Retrieved 2026-02-25, License - CC BY-SA 3.0
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}


cv_plots <- function(result, f_true, cutoff=1e3, lower=1e-16, upper=1e6){


  meanplt <- result$vh_mean %>% as.data.frame() %>% filter(lambda>= lower & lambda <= upper) %>% ggplot(aes(x=lambda, y=vh_mean)) + geom_point()+ scale_x_log10() + theme_bw()+
    theme(panel.grid=element_blank()) + ylab("Volume-averaged conditional squared residual") + xlab("λ")+ ggtitle("(ii)") +
    geom_vline(xintercept=cutoff, col="black",linetype ="dashed") + theme(aspect.ratio=1)

  sdplt <- result$vh_sd %>% as.data.frame() %>%  filter(lambda>= lower & lambda <= upper) %>% ggplot(aes(x=lambda, y=vh_sd)) + geom_point()+ scale_x_log10() + theme_bw()+
    theme(panel.grid=element_blank())+ ylab("Standard Deviation") + xlab("λ")+ ggtitle("(iii)") + theme(aspect.ratio=1)


  S <- seq(0,0.999,length.out=500)
  # Setup Data generating curve
  #f_true <- flower1d2D(500,5,noise=0, t=S)$X
  df <- data.frame(t=S, x=f_true[,1],y=f_true[,2], repetition = 0)

  lams <- result$cv[,1]
  lams_sub <- lams[lams < cutoff]
  cv_vals <- result$cv[,2]
  cv_sub <- cv_vals[lams < cutoff]
  #localmins <- find_peaks(-cv_vals)
  lam_optim <- lams[which.min(cv_vals)]
  lam_optim_sub <- lams_sub[which.min(cv_sub)]

  optimal_ind <- min(which(lams == lam_optim_sub))


  #localmins_cut <- localmins[localmins %in% which(lams< 1e-4 & lams> 1e-9)]
  #df_cut <- data.frame(candidate=1:length(localmins_cut), lambda = lams[localmins_cut])
  cvplt <- result$cv %>% as.data.frame() %>% filter(lambda>= lower & lambda <= upper) %>% ggplot(aes(x=lambda, y=cv))  +  geom_point()+ scale_x_log10() + theme_bw()+
    theme(panel.grid=element_blank()) + ylab("Coefficient of Variation") + xlab("λ")+ ggtitle("(i)") +
    geom_vline(aes(xintercept=lam_optim_sub),color="blue" ,lty="dashed") + theme(aspect.ratio=1)

  spl_optim_eval <- result$spline_list[[optimal_ind]](S)
  df_fit_lamb <- data.frame(t = S, x = spl_optim_eval[,1], y = spl_optim_eval[,2], group ="Estimate")
  df_fit_lamb <- rbind(df_fit_lamb, data.frame(t = S, x = f_true[,1], y = f_true[,2], group = "True"))
  df_fit_lamb <- rbind(df_fit_lamb, data.frame(t = -1, x=result$dataset[,1], y=result$dataset[,2], group="Data"))
  df_fit_lamb <- df_fit_lamb %>% mutate(group = factor(group, levels=c("Data","True","Estimate")))

  #df_fit_lamb <- do.call(rbind, lapply(seq_along(localmins_cut), function(i){
  #  t = S
  # output <- result$spline_list[[localmins_cut[i]]](S)
  # lambda <- lams[localmins_cut[i]]

  # as.data.frame(cbind(t = t, x = output[,1], y = output[,2], lambda = lambda))
  #})) %>% mutate(group = "Sample")

  labels <- c("Data", "Estimate", "True")
  (plotFit <- df_fit_lamb %>%  ggplot()+
      # Add data points
      geom_point(data=subset(df_fit_lamb, group=="Data"),aes(x=x,y=y, col=group, alpha=group, linetype=group),size=0.7)+
      geom_path(data = subset(df_fit_lamb, group != "Data"),aes(x=x,y=y, col = group, linetype = group, alpha=group),linewidth=1)+

      scale_alpha_manual(guide="legend", name ="group", breaks=labels, values = c(0.2, 0.8, 0.8))+
      scale_color_manual(guide="legend", name ="group", breaks =labels, values = c("black","blue","red") )+
      scale_linetype_manual(guide="legend",name ="group", breaks = labels,values =c("blank", "dashed", "solid"))+
      labs(col="",linetype="", alpha="")+
      theme_bw()+
      theme(panel.grid=element_blank())+
      theme(legend.background=element_rect(fill="white",color="black"),
            legend.title=element_blank()) + ggtitle("(iv)") + theme(legend.position="bottom", aspect.ratio=1)
  )

  cvplt + meanplt + sdplt + plotFit + plot_layout(ncol=4)
}

# periodic
f_true <- flower1d2D(1000,5,noise=0, t=seq(0,1,length.out=1000))$X
(plt_per<-cv_plots(pme_periodic, f_true, cutoff=1e-4))
ggsave(file.path(figure_dir, "cv_periodic_example.png"),plt_per, width = 7200, height = 2400, units="px",dpi=600)

# interval
f_true <- halfflower1d2D(1000, petals=0, rotation= 3*pi/2, noise=0.0)$X

(plt_int<-cv_plots(pme_interval, f_true, cutoff=1e-4))
ggsave(file.path(figure_dir, "cv_interval_example.png"),plt_int,width = 7200, height = 2400, units="px",dpi=600)




