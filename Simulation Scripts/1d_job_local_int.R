library(vegan) # isomap
library(plotly) # plotting
library(foreach) # parallelization
library(doParallel)
library(Riemann) # for 2d variance heterogeneity
library(pracma)
library(TDA)
library(splines)

#library(rgl)
#dos2unix
#squeue --start --me

#==== HELPER FUNCTIONS ==================
#cooltools functions necessary. Taken from github to avoid loading unnecessary functions
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
spline_error <- function(Xdata, t, splineresult){
  mean((rowSums( (Xdata-splineresult(t))^2)))
}

# Helper function for variance heterogeneity given smoothed residuals from local_var_1d
var_het <- function(sigma2_hat) {
  mean_sigma2 <- mean(sigma2_hat)
  sd_sigma2   <- sd(sigma2_hat)
  cv          <- sd_sigma2 / (mean_sigma2 + 1e-12)  # coefficient of variation
  spread      <- max(sigma2_hat) - min(sigma2_hat)
  list(mean = mean_sigma2, sd = sd_sigma2, cv = cv, spread = spread)
}


#  Compute smoothed residuals using kernel smoothing
local_var_1d <- function(t, resid_sq, periodic=F,
                         grid_eval = 100,
                         bandwidth = NULL,
                         kernel = "normal") {
  ## periodic extension: t-1, t, t+1
  t_ext <- t %% 1 # CHANGE MENG added %% 1
  e2_ext <- resid_sq
  
  # If periodic (set periodic=T for circle template)
  if(periodic){
    t_ext  <- c(t_ext - 1, t_ext, t_ext + 1)
    e2_ext <-  rep(resid_sq, 3)
    print("Doing periodic")
  }
  
  
  
  
  ## ---------- automatic bandwidth selection ----------
  if (is.null(bandwidth)) {
    bw <- NA_real_
    
    if (requireNamespace("KernSmooth", quietly = TRUE)) {
      ord <- order(t)
      tx  <- t[ord]
      ey  <- resid_sq[ord]
      bw_try <- try(KernSmooth::dpill(x = tx, y = ey), silent = TRUE)
      
      if (!inherits(bw_try, "try-error") && is.finite(bw_try) && bw_try > 0) {
        bw <- bw_try
      }
    }
    
    if (!is.finite(bw) || bw <= 0) {
      bw <- stats::bw.SJ(t)
    }
    
    bandwidth <- bw
  }
  ## ---------------------------------------------------
  # grid <- seq(0, 1, length.out = grid_eval) CHANGE MENG
  grid <- seq(0, 1, length.out = grid_eval + 1)[-(grid_eval + 1)]  # [0,1)
  
  ks <- ksmooth(
    x        = t_ext,
    y        = e2_ext,
    kernel   = kernel,
    bandwidth = bandwidth,
    x.points = grid
  )
  sigma2_hat <- pmax(ks$y, 0)  # clamp to nonnegative
  list(s = grid, sigma2_hat = sigma2_hat, bandwidth = bandwidth)
}

# Circular projection, for initialization puprposes
project_circle <- function(X){
  Xc <- scale(X, scale=F)
  (atan2(Xc[,2],Xc[,1])/(2*pi)) %% 1
}

# Optimization over a grid. This is to find the best projection index on a 1d Manifold
project_grid <- function(Xraw, res, gridSize = 1000) {
  # Initialize grid
  t_grid <- seq(0, 1, length.out = gridSize)
  d <- ncol(Xraw)
  
  # Define curve based on spline result and vectorize to find min
  curve <- res(t_grid)
  idx <- apply(Xraw, 1, function(x) {
    which.min(rowSums((curve - matrix(rep(x, each=gridSize), nrow=gridSize, ncol = d))^2))
  })
  
  # Return ordering
  t_grid <- t_grid[idx]
  #cbind(Xraw[t_sort$ix, ],t_sort$x)
}

# Find projection index using an optimization routine. Grid method is preferred for global solutions
project_optimize <- function(Xraw, res, periodic = F) {
  # Intialize t
  n <- nrow(X)
  topt <- rep(0, n)
  # Apply iterative optimization scheme based on spline fit
  
  for (i in 1:n) {
    ff <- function(t) {
      sum((Xraw[i, ] - res(t) )^2)
    }
    t[i] <- optimize(ff, interval = c(0, 1))$minimum
  }
  # Return t
  ret <- topt
  
  if(periodic)
    ret <- topt %% 1
  
  return(ret)
}

# An isomap method to find the smallest k so that the graph is not fragmented
iso_try <- function(X, k=1){
  k<-k
  repeat{
    output <- try(isomap(dist(X),ndim=1,k=k),silent = T)
    if(!inherits(output, "try-error")) return(list(iso=output, k = k))
    
    k <- k+1
    
  }
}



# SPLINES AND PME: INTERVAL ============================
normalize <- function(t){
  a <- min(t)
  b <- max(t)
  (t-a)/(b-a)
}

# The reproducing kernel for the interval
Rk <- function(s,t){
  1/3*pmin(s,t)^3 - 1/2 *pmin(s,t)^2 * (s+t) + pmin(s,t)*s*t
}

# The exact spline on the interval
spline1d_interval <- function(X,t,lambda){
  n <- nrow(X)
  # T matrix
  t <- as.vector(t)
  Tmat <- cbind(1, t)
  
  #Sigma <- combn(1:n, n, Vectorize(function(i, j) Rk(t[i],t[j])))[,,1]
  Sigma <- outer(t, t, function(s, t) Rk(s, t))
  
  M <- Sigma + n*lambda*eye(n)
  
  # Theta factor (multiplies y) same for x and y coordinate. Calculate once.
  #Minv <- solve(M)
  # K is symmetric positive definite by propertie of rk. M is also symmetric positive definite
  # chol2inv chol M is 2x faster than normal solve
  Minv <- chol2inv(chol(M))
  #Minv <- solve(M)
  
  # Precompute coefficients
  theta_coeffs <- ( solve(t(Tmat) %*% Minv %*% Tmat) %*% t(Tmat) %*% Minv) %*% X # coefficients for null space
  alpha_coeffs <- Minv %*% (X - Tmat %*% theta_coeffs)  # (n x D) # coefficients for penalized space
  
  # calculate periodic function each dimension. Help of chatgpt ===
  fx <- function(x) {
    Fmat <- outer(x, t, function(xi, tj) Rk(xi, tj))  # length(x) x n
    Fmat %*% alpha_coeffs + cbind(1, x) %*% theta_coeffs
  }
  
  # is this correct?
  attr(fx, "alpha") <- sum(diag(t(alpha_coeffs) %*% Sigma %*% alpha_coeffs))
  
  return(fx)
}

# The approximate spline on the interval
spline1d_interval_approx <- function(X,t,lambda){
  
  # Get spline fits for arbitrary ambient dimension
  D <- ncol(X)
  N <- nrow(X)
  fits <- vector(mode="list", length=D)
  
  fits <- apply(X,2, function(data_col){
    smooth.spline(t, data_col,lambda=lambda)
  })
  
  
  # Calculate penalty
  # Return penalty
  penalty <- 0
  for(i in seq_along(fits)){
    fiti <- fits[[i]]
    c <- fiti$fit$coef
    knotvals <- fiti$fit$knot
    
    t_grid <- seq(min(knotvals), max(knotvals), length.out=10*N)  # fine grid
    delta <- diff(t_grid)[1]  # step size for numerical integration
    
    B_dd <- splineDesign(
      knots = knotvals,
      x = t_grid,
      ord = 4,        # cubic spline
      derivs = 2,     # second derivatives
      outer.ok = TRUE
    )
    
    Sigma <- t(B_dd) %*% B_dd * delta
    
    penalty <- penalty + t(c) %*% Sigma %*% c
  }
  
  
  spline_fun <- function(tnew) {
    eval <- lapply(fits, function(spline) predict(spline, tnew)$y)
    do.call(cbind, eval)
  }
  
  
  
  
  attr(spline_fun, "alpha") <- penalty 
  spline_fun
  
}


# ========== Error for fit ====================
spline_error <- function(Xdata, t, splineresult){
  mean((rowSums( (Xdata-splineresult(t))^2)))
}

# PME method
# Parameters mostly identical to above.
pmeClosed_1d_interval <- function(data, lambdas, gridSize = 1000, max_iter = 10, optimize_lambda=F,
                                  err_tol=0.01, plot_lambda=T,lvp_grid=NULL, cl=NULL,parallel=F, k_init = 1,
                                  approx=F, min_iter=1)
{
  n <- nrow(data)
  X <- data
  D <- ncol(data)
  
  lvp_grid <- ifelse(is.null(lvp_grid), 10*n, lvp_grid)
  
  # Initialize step: project onto unit circle. This is independent of the lambda choice.
  tinit <- NULL
  
  # Rule of thumb for choosing iso. Smallest possible without error
  isores <- iso_try(X, k_init)
  iso <- isores$iso
  k <- isores$k
  print(paste0("k is ", toString(k)))
  
  tinit <- normalize(iso$points)
  
  # Projection wrapper
  miniproj <- function(Xraw, newSpline){
    project_grid(Xraw, newSpline, gridSize)
  }
  
  spline1d_interval_func <- ifelse(approx, spline1d_interval_approx,
                                   spline1d_interval)
  
  # Setup parallel cluster
  # Setup cluster if not provided
  stop_after <- FALSE
  if (is.null(cl) && parallel) {
    numCores <- parallel::detectCores() - 1
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    stop_after <- TRUE
  }
  # run parallel? T or F
  operator <- if(parallel) `%dopar%` else `%do%` #thanks chatgpt
  
  
  results <-operator( foreach(i = 1:length(lambdas), 
                              .packages = c("stats", "pracma", "splines"),
                              .export = c("spline1d_interval","project_grid", "project_optimize",
                                          "spline_error","Rk", "local_var_1d", "var_het","normalize", "splineDesign")), {
                                            lambda <- lambdas[i]
                                            #splines <- vector(mode="list",length=max_iter)
                                            alphas <-rep(0,max_iter)
                                            fidelities <- rep(0, max_iter)
                                            
                                            # Setup first spline fit
                                            spl <- spline1d_interval_func(X,tinit, lambdas[i])
                                            alphas[1] <- sum(attr(spl, "alpha"))
                                            
                                            
                                            # Setup first projection 
                                            # tt <- list() # list of t values
                                            t_proj<- normalize(miniproj(X,spl)) # evaluate first projection set
                                            fidelities[1] <- spline_error(X, t_proj, spl)
                                            
                                            # Repeat until error small
                                            err <- err_tol + 1
                                            iter <- 1
                                            #Diter <- rep(Inf, max_iter)
                                            # Initial error
                                            #Diter[iter] <- spline_error(X, t_proj, spl)
                                            
                                            
                                            # Iteriate until error convergences
                                            while( (err > err_tol || iter < min_iter) &&  iter < max_iter){
                                              iter <- iter+1
                                              
                                              #  Get newest spline fit
                                              spl <- spline1d_interval_func(X, t_proj, lambdas[i])
                                              alphas[iter] <- attr(spl, "alpha")
                                              
                                              # Calculate new projection
                                              t_proj <- normalize(miniproj(X,spl)) # evaluate first projection set
                                              fidelities[iter] <- spline_error(X,t_proj, spl)
                                              
                                              # Get newest error
                                              #Diter[iter] <- spline_error(X,tt[[iter]], splines[[iter]])
                                              
                                              # Compare old error
                                              err <-abs((fidelities[iter] - fidelities[iter-1])/fidelities[iter-1])
                                            }
                                            
                                            resid_sq <- rowSums( (X - spl(t_proj))^2)
                                            
                                            # Compute variance heterogeneity
                                            vh <- NULL
                                            if(optimize_lambda){
                                              lvp <- local_var_1d(t_proj, resid_sq,grid_eval=lvp_grid)$sigma2_hat
                                              vh <- var_het(lvp) # get coefficient of variation
                                            }
                                            
                                            
                                            list(spline = spl, X = X, tproj = t_proj, lambda = lambdas[i],
                                                 resid_sq = resid_sq, cv = vh$cv,alphas=alphas[1:iter],
                                                 mean_sq = fidelities[1:iter], cost_functional = lambdas[i]*alphas[1:iter] + fidelities[1:iter], vh_mean = vh$mean,
                                                 vh_sd = vh$sd)
                                          })
  # Aggregate results and find minimizing spline
  spline_list <- lapply(results, function(res) res$spline)
  lambda_list <- lambdas
  alpha_list <- lapply(results, function(res) res$alphas)
  mean_sq_list <- lapply(results, function(res) res$mean_sq)
  cost_functional_list <- lapply(results, function(res) res$cost_functional)
  
  
  list_return <- list(spline_list = spline_list, lambda_list = lambda_list,
                      alpha_list = alpha_list, mean_sq_list = mean_sq_list,
                      cost_functional_list = cost_functional_list, dataset=data)
  
  
  # If optimize lambda true, add optimal spline by cv.
  # also add cv list for reference
  if(optimize_lambda){
    cv_list <- sapply(results, function(res) res$cv)
    vh_mean_list <- sapply(results, function(res) res$vh_mean)
    vh_sd_list <- sapply(results, function(res) res$vh_sd)
    
    min_ind <- which.min(cv_list)
    list_return$spline_optimal = spline_list[[min_ind]]
    list_return$lambda_optimal = lambda_list[min_ind]
    list_return$cv = cbind(lambda = lambdas, cv = cv_list)
    list_return$vh_mean = cbind(lambda = lambdas, vh_mean = vh_mean_list)
    list_return$vh_sd = cbind(lambda = lambdas, vh_sd = vh_sd_list)
  }
  # Return list structure
  #list_return <- list(spline_optimal = spline_list[[min_ind]], lambda_optimal = lambda_list[min_ind],
  #  spline_list = spline_list, lambda_list = lambda_list, cv_list = cv_list)
  
  
  # stop the cluster
  on.exit({
    if (stop_after) stopCluster(cl)
  }, add = TRUE)
  
  return(list_return)
}



#========== datasets ============


halfflower1d2D <- function(n, petals, noise, t, rotation=0){
  t0 <- if(missing(t)) sort(runif(n)) else t
  
  r <- (1 + 0.3 * sin(petals *pi* t0))
  return(list(t=t0, X = cbind(r * cos(pi*t0),  r * sin(pi*t0)) %*% rotation2(rotation) +rnorm(2*n, sd=noise)))
} 

#======tda========
runTDA_1d2D <- function(data, alpha=0.01, boot_size=30, by = 0.01)
{
  Xlim <- c(min(data[,1]), max(data[,1]))
  Ylim <- c(min(data[,2]), max(data[,2]))
  lim <- cbind(Xlim, Ylim)
  
  Diag <- gridDiag(data, distFct, lim = lim, by = by)
  boot <- bootstrapDiagram(X=data, FUN=distFct, lim = lim,by=by, parallel=T,alpha=alpha,B=boot_size)
  
  list(Diag=Diag, boot=boot, alpha =alpha)
}




# ============================ RUN SIMULATION =================================
# Deliverables:
# A: TDA plot + fits + functional
# B: CV plot 
# 1d curve in 2D



# =============== 1d ======================
# 1. load data
set.seed(1)
N <-2500
boot <- 100
lambda_N <- 150
max_it <-50
min_it <- 10


#dataset_interval <- spiral1d2D_arclength(N, period=5*pi, noise=0.02)
dataset_interval <- halfflower1d2D(N, petals=0, rotation= 3*pi/2, noise=0.02)

# 2. perform tda and save result
tda_interval <- runTDA_1d2D(dataset_interval$X, boot_size=boot)
saveRDS(tda_interval, "results_1d/1d_interval_tda_local_exact.rds")

# 3. perform pme and save result
lambdas_int <- 10^(seq(-12, 3, length.out=lambda_N))
pme_interval <- pmeClosed_1d_interval(dataset_interval$X, lambdas_int, err_tol = 0.001, max_iter = max_it,  min_iter = min_it, optimize_lambda=T,
                                      parallel=T, gridSize=N,approx=F, lvp_grid=10*N)

saveRDS(pme_interval, "results_1d/1d_interval_pme_local_exact.rds")