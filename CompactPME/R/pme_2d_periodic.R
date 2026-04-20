# Generated from create-CompactPME.Rmd: do not edit by hand

#' # The spline fit for the 1d Interval
#' 
#' @param data an Nx3 data matrix where N is the number of observations
#' @param lambdas A numeric vector of lambda values where we will perform a spline fit on. Lambda values typically are in the range of 1e-12 to 1e-2 for interpolation and for limiting lambda form respectively.
#' @param min_iter the minimum number of total iterations of PME. Default is 1 for speed, but choosing a higher number is useful if you want to visualize the cost function.
#' @param max_iter the maximum number of iterations of PME. Prevents the algorithm from stalling.
#' @param err_tol error tolerance for convergence of the algorithm. Lower values result in longer runtimes, and more likely to hit the maximum iteration cap.
#' @param optimize_lambda whether or not to perform the method of section 5, Tuning Parameter Selection. If true, additional values will be returned, being the coefficient of variation at each lambda, along with the mean residual and standard deviation residual used to calculate it.
#' @param lvp_grid The number denoting the size of the grid used for kernel smoothing for the tuning parameter selection. If left as NULL, it will choose a grid of 10 times the sample size.
#' #' @param cl A a makeCluster object for parallelization. For most cases, should be left null so that a new cluster is made upon runtime.
#' @param parallel whether or not to run the algorithm parallel. The algorithm will be parallelized with respect to lambda and is recommended if your lambdas vector contains many values. Not recommended for only a handful of lamda values unless your dataset is large.
#' @param k_init The starting value for a search of the ISOMAP nearest neighbor parameter k. Our method is to find the smallest such k such that the ISOMAP graph is connected, as detailed in the iso_try function. The default is 1 to find the smallest possible.
#' @param approx whether or not to use the approximate spline method. We recommend to set this to be TRUE for very large lambda vectors and or very large sample sizes.
#' @export
#' @importFrom foreach %dopar% %do%
pme_2d_periodic <- function(data, lambdas, min_iter=1, max_iter = 10, optimize_lambda=F, err_tol = 0.01,parallel=F, cl = NULL, init_iso=F,k_init =10,
                           resid_smooth_size = 1000,bw_grid = 1000, pelletier_bw=0.1) {
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
                                 .packages = c("stats", "pracma", "cooltools","Riemann"),
                                 .export = c("spline2d","project_optimize2",
                                             "spline_error","qm","expr","Rpp","m", "local_var_2d", "var_het",
                                             "d_circ_vec","pelletier_kernel")), {
                                               
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
                                                 vhi <- var_het(lvp$s2_intrin) # get coefficient of variation
                                                 vhe <- var_het(lvp$s2_extrin) 
                                                 vhp <- var_het(lv_pelletier)
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
