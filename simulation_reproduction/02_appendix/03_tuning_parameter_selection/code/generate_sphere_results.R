library(vegan) # isomap
library(plotly) # plotting
library(foreach) # parallelization
library(doParallel)
library(Riemann) # for 2d variance heterogeneity
library(pracma)
library(TDA)
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



iso_try <- function(X, k=1, d=1){
  k<-k
  repeat{
    output <- try(vegan::isomap(dist(X),ndim=d,k=k),silent = T)
    if(!inherits(output, "try-error")) return(list(iso=output, k = k))
    
    k <- k+1
    
  }
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

# =============== 2d Manifolds =============== 

# Some variance heterogeneity methods===============
var_het <- function(sigma2_hat,eps_adj=0) {
  mean_sigma2 <- mean(sigma2_hat)
  sd_sigma2   <- sd(sigma2_hat)
  cv          <- sd_sigma2 / (mean_sigma2 + eps_adj)  # coefficient of variation
  spread      <- max(sigma2_hat) - min(sigma2_hat)
  list(mean = mean_sigma2, sd = sd_sigma2, cv = cv, spread = spread)
}


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
d_circ_vec <- function(x_new, x_data) {
  # x_new: vector (length d)
  # x_data: matrix n_train x d
  dots <- as.vector(x_data %*% x_new)  # fast matrix multiply in R
  pmax(pmin(dots,1),-1) |> acos()
}

#' Pelletier, B. (2006). Non-parametric regression estimation on closed Riemannian manifolds. Journal of Nonparametric Statistics, 18(1), 57-67.
pelletier_kernel <- function(x_news, x_data, y_data, h) {
  sapply(1:nrow(x_news), function(i){
    rho <- d_circ_vec(x_news[i,], x_data)
    k <- dnorm(rho / h)*ifelse(rho<1e-16, 1, rho/sin(rho))
    sum(k * y_data) / sum(k)
  })
}


# Splines and intialization functions==================
# Kernel functions
m <- 2

# qm function
qm <- function(z){
  # if z>=1, issue
  # gpt suggestion. z forced to be between -1 and 1e-12
  z <- pmin(pmax(z,-1), 1-1e-12)
  
  
  zinv <- 1-z
  0.5*(log(1+sqrt(2/zinv))*( 3*zinv^2 - 2*zinv)-
         12*(zinv/2)^(1.5) + 6*(zinv/2)+1)
}
# kernelhelper
expr <- function(z) { 1/(4*pi) *(qm(z) - 1/3)}


# Kernel function
Rpp <- function(v1,v2){
  expr(sum(v1*v2))
} 

# error helper
spline_error <- function(Xdata, t, splineresult){
  mean((rowSums( (Xdata-splineresult(t))^2)))
}

# Spline fit helper
spline2d <- function(data,Tt, lambda){
  n <- nrow(data)
  #Form R kernel matrix (chatgpt)
  # Normalize rows of data
  #Znorm <- data / sqrt(rowSums(data^2))  # n x 3
  #Tt <- scale(Tt,scale=F)
  Znorm <- Tt / sqrt(rowSums(Tt^2))
  # Compute all pairwise cos(gamma) = dot products
  cosGamma <- Znorm %*% t(Znorm)  # n x n
  
  #View(cosGamma)
  # Apply kernel function to all pairs at once
  Rmat2 <- matrix(0, nrow=n, ncol=n)
  # Diagonal: limit value
  diag(Rmat2) <- 1/(24*pi)
  
  # Off-diagonal: apply expr to all non-diagonal entries
  diag(cosGamma) <- NA  # ignore diagonal
  Rmat2[!is.na(cosGamma)] <- expr(cosGamma[!is.na(cosGamma)])
  
  #============== Precompute matrices
  m <- 2
  r <- 0:m
  Tmat <- ones(n,1)
  M <- Rmat2 + n * lambda * eye(n)
  Minv <- solve(M)
  
  # constant part
  const_factor <- solve( t(Tmat) %*% Minv %*% Tmat) %*% t(Tmat) %*% Minv %*% data
  coeffs <- Minv %*% (data - Tmat %*% const_factor) 
  
  
  #================== Return spline function
  f_eval_fast <- function(Pnew) {
    Pnew <- matrix(Pnew,ncol=3)
    
    # normalize rows
    Pnew_norm <- Pnew / sqrt(rowSums(Pnew^2))
    Z_norm <- Tt / sqrt(rowSums(Tt^2))
    
    # compute all pairwise dot products (m x n)
    cosGamma <- Pnew_norm %*% t(Z_norm)
    
    # apply kernel function elementwise
    # cosGamms is all dot products. We evaluate R kernel separately whether
    # cosGamma is 1 or not
    Rnew <- cosGamma
    oneloc <- which(cosGamma>=1, arr.ind=T)
    cloc <- which(cosGamma < 1, arr.ind=T)
    Rnew[oneloc] <- 1/(24*pi)
    Rnew[cloc] <- expr(Rnew[cloc])
    #new
    
    #Rnew <- expr(cosGamma)
    # evaluate spline
    Rnew %*% coeffs + matrix(rep(const_factor, each = nrow(Pnew)), nrow = nrow(Pnew))
  }
  
  # add penalty attribute
  attr(f_eval_fast, "alpha") <-sum(diag(t(coeffs) %*% Rmat2 %*% coeffs))
  # Return function
  f_eval_fast
  
  #return(Rmat2)
}

# Projection via optimization
# Optimize surface based on spline values and previous projection parameters
# (cartesian)
project_optimize2 <- function(X, f, init) {
  param_return <- matrix(0, nrow = nrow(X), ncol = 2)
  
  init_sph <- cart2sph(init)[,-3]
  
  objective <- function(param,Xi) {
    input <- matrix(param, ncol=2)
    input_cart <- sph2cart(cbind(input, 1))
    sum((Xi - f(input_cart)  )^2  )
  }
  
  for (i in 1:nrow(X)) {
    nl <- nlm(objective, init_sph[i,], Xi=X[i,])$estimate
    param_return[i,] <- nl
  }
  sph2cart(cbind(param_return, 1))
}

# PME for the 2-sphere =================
pme_2d_periodic <- function(data, lambdas, min_iter=1, max_iter = 10, optimize_lambda=F, err_tol = 0.01,parallel=F, cl = NULL, init_iso=F,k_init =10,
                            resid_smooth_size = 1000,bw_grid = 1000, pelletier_bw=0.1, eps_adj=0) {
  # Setup data and perform isomap to get initial parameterization
  X <- data
  m <- 2
  if(init_iso){
    isores <- iso_try(X, k_init,3)
    iso3 <- isores$iso
    print(paste0("k is ", toString(isores$k)))
    
    t3 <- scale(iso3$points,scale=F)
    Znorm <-t3 / sqrt(rowSums(t3^2))
  }else{
    datac <- scale(X,scale=F)
    Znorm <- datac / sqrt(rowSums(datac^2))
  }
  # Setup cluster if not provided
  
  stop_after <- FALSE
  if (is.null(cl) && parallel) {
    numCores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(numCores)
    doParallel::registerDoParallel(cl)
    stop_after <- TRUE
  }
  # run parallel? T or F
  operator <- if(parallel) `%dopar%` else `%do%`
  
  
  # Parallel loop over lambdas
  results <- operator( foreach::foreach(i = 1:length(lambdas), 
                                        .packages = c("stats", "pracma","Riemann"),
                                        .export = c("spline2d","project_optimize2",
                                                    "spline_error","qm","expr","Rpp","m", "local_var_2d", "var_het",
                                                    "d_circ_vec","pelletier_kernel", "fibonaccisphere")), {
                                                      
                                                      Xi <- X
                                                      Znorm_i <- Znorm
                                                      
                                                      lambda <- lambdas[i]
                                                      
                                                      # First spline
                                                      #splines <- list()
                                                      alphas <-rep(0,max_iter)
                                                      spl <- spline2d(Xi, Znorm_i, lambda)
                                                      alphas[1] <- attr(spl, "alpha")
                                                      projection <- project_optimize2(Xi, spl, Znorm_i)
                                                      
                                                      # Iterative refinement
                                                      err <- err_tol + 1
                                                      iter <- 1
                                                      fidelities <- rep(Inf, max_iter)
                                                      fidelities[iter] <- spline_error(Xi, projection, spl)
                                                      
                                                      #t1 <- proc.time()
                                                      while( (err > err_tol || iter < min_iter) &&  iter < max_iter){
                                                        iter <- iter + 1
                                                        spl <- spline2d(Xi, projection, lambda)
                                                        alphas[iter] <- attr(spl, "alpha")
                                                        projection <- project_optimize2(Xi, spl, projection)
                                                        fidelities[iter] <- sum( rowSums((Xi - spl(projection))^2))
                                                        err <- abs((fidelities[iter] - fidelities[iter - 1]) / fidelities[iter - 1])
                                                      }
                                                      
                                                      #print(paste0("Lambda ", toString(lambda), ": iteration time: ",toString(proc.time()-t1)))
                                                      
                                                      
                                                      resid_sq <- rowSums( (X - spl(projection))^2)
                                                      
                                                      # Compute variance heterogeneity, intrinsic and extrinsic
                                                      #t1 <- proc.time()
                                                      lvp <- NULL
                                                      vhi <- NULL
                                                      vhe <- NULL
                                                      vhp <-NULL
                                                      
                                                      # If optimizing lambda, perform variance heterogeneity routine
                                                      if(optimize_lambda){
                                                        lvp <- local_var_2d(projection, resid_sq, eval_size=resid_smooth_size, bw_grid=bw_grid)
                                                        lv_pelletier <- pelletier_kernel(fibonaccisphere(resid_smooth_size),
                                                                                         projection, resid_sq, pelletier_bw)
                                                        #print(lv_pelletier)
                                                        # Get coefficient of variation
                                                        vhi <- var_het(lvp$s2_intrin,eps_adj) # get coefficient of variation
                                                        vhe <- var_het(lvp$s2_extrin,eps_adj) 
                                                        vhp <- var_het(lv_pelletier,eps_adj)
                                                      }
                                                      
                                                      
                                                      list(spline = spl, tproj = projection, lambda = lambdas[i],
                                                           resid_sq = resid_sq, cv_intrinsic = vhi$cv, cv_extrinsic = vhe$cv, cv_pelletier = vhp$cv, 
                                                           mean_intrinsic = vhi$mean, mean_extrinsic = vhe$mean, mean_pelletier = vhp$mean,
                                                           sd_intrinsic = vhi$sd, sd_extrinsic = vhe$sd, sd_pelletier = vhp$sd,
                                                           bw_intrinsic =  lvp$bw_intrin, bw_extrinsic = lvp$bw_extrin, bw_pelletier = pelletier_bw,
                                                           alphas=alphas[1:iter], mean_sq = fidelities[1:iter],
                                                           cost_functional = lambdas[i]*alphas[1:iter] + fidelities[1:iter])
                                                    })
  
  
  
  # Aggregate results and find minimizing spline
  spline_list <- lapply(results, function(res) res$spline)
  lambda_list <- lambdas
  proj_list <- lapply(results, function(res) res$tproj)
  resid_sq_list <- lapply(results, function(res) res$resid_sq)
  cost_functional_list <- lapply(results, function(res) res$cost_functional)
  
  alpha_list <- lapply(results, function(res) res$alphas)
  mean_sq_list <- lapply(results, function(res) res$mean_sq)
  
  
  list_return <- list(spline_list = spline_list, proj_list = proj_list, resid_sq_list = resid_sq_list, lambda_list = lambda_list,
                      alpha_list = alpha_list, mean_sq_list = mean_sq_list,
                      cost_functional_list = cost_functional_list,dataset=data)
  
  # If we optimize lambda, add optimal spline data
  if(optimize_lambda){
    
    # intrinsic
    cvi_list   <- sapply(results, function(res) res$cv_intrinsic)
    meani_list <- sapply(results, function(res) res$mean_intrinsic)
    sdi_list   <- sapply(results, function(res) res$sd_intrinsic)
    
    # extrinsic
    cve_list   <- sapply(results, function(res) res$cv_extrinsic)
    meane_list <- sapply(results, function(res) res$mean_extrinsic)
    sde_list   <- sapply(results, function(res) res$sd_extrinsic)
    
    # pelletier
    cvp_list   <- sapply(results, function(res) res$cv_pelletier)
    meanp_list <- sapply(results, function(res) res$mean_pelletier)
    sdp_list   <- sapply(results, function(res) res$sd_pelletier)
    
    
    # which spline minimizes CV?
    mini_ind <- which.min(cvi_list)
    mine_ind <- which.min(cve_list)
    minp_ind <- which.min(cvp_list)
    
    
    # store optimal splines
    list_return$spline_optimal_intrinsic  <- spline_list[[mini_ind]]
    list_return$lambda_optimal_intrinsic  <- lambda_list[mini_ind]
    
    list_return$spline_optimal_extrinsic  <- spline_list[[mine_ind]]
    list_return$lambda_optimal_extrinsic  <- lambda_list[mine_ind]
    
    list_return$spline_optimal_pelletier  <- spline_list[[minp_ind]]
    list_return$lambda_optimal_pelletier  <- lambda_list[minp_ind]
    
    list_return$cv <- cbind(
      lambda     = lambda_list,
      intrinsic  = cvi_list,
      extrinsic  = cve_list,
      pelletier  = cvp_list
    )
    
    list_return$mean <- cbind(
      lambda     = lambda_list,
      intrinsic  = meani_list,
      extrinsic  = meane_list,
      pelletier  = meanp_list
    )
    
    list_return$sd <- cbind(
      lambda     = lambda_list,
      intrinsic  = sdi_list,
      extrinsic  = sde_list,
      pelletier  = sdp_list
    )
    
  }
  
  # stop the cluster
  on.exit({
    if (stop_after) parallel::stopCluster(cl)
  }, add = TRUE)
  
  list_return$template <- "S2"
  list_return$D <- D
  class(list_return) <- "cpme"
  
  return(list_return)
}

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
  boot <- bootstrapDiagram(X=data, FUN=distFct, lim = lim,by=by, parallel=T,alpha=alpha,B=boot_size)
  
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
  boot <- bootstrapDiagram(X=data, FUN=distFct, lim = lim,by=by, parallel=T,alpha=alpha,B=boot_size)
  
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
N <-250
dataset <- flower2d3D_unif(N, petals=5, noise=0.02)
#plot_data_3D(dataset$X)

# 2. TDA step
#boot <- 100
#tda_sphere <- runTDA_2d3D(dataset$X, boot_size=boot)
#saveRDS()
#tda_2_plot_3D(tdares)
#tdapath <- paste0("results_sph/flower_tda", job_label, ".rds")
#saveRDS(tda_sphere, tdapath)

# 3. Run PME
lambda_N <-10
lambdas <- 10^(seq(-14,-3,length.out=lambda_N))
pme_sphere <- pme_2d_periodic(dataset$X, lambdas, optimize_lambda = T,parallel=T, resid_smooth_size = 5*N)
pmepath <- paste0("results_sph/flower_pme", job_label, ".rds")
saveRDS(pme_sphere, pmepath)


N <-250
dataset <- unifmoon2d3D(N,noise=0.02)
#plot_data_3D(dataset$X)

# 2. TDA step
boot <- 10
tda_sphere <- runTDA_2d3D(dataset$X, boot_size=boot)
#saveRDS()
#tda_2_plot_3D(tdares)
tdapath <- paste0("results_sph/cashew_tda", job_label, ".rds")
saveRDS(tda_sphere, tdapath)

# 3. Run PME
lambda_N <-10
lambdas <- 10^(seq(-14,-3,length.out=lambda_N))
pme_sphere <- pme_2d_periodic(dataset$X, lambdas, optimize_lambda = T,parallel=T,resid_smooth_size = 5*N)
pmepath <- paste0("results_sph/cashew_pme", job_label, ".rds")
saveRDS(pme_sphere, pmepath)

