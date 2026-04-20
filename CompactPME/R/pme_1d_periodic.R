# Generated from create-CompactPME.Rmd: do not edit by hand

#' # The spline fit for the 1d Interval
#' 
#' @param data an NxD data matrix where N is the number of observations and D is either 2 or 3.
#' @param lambdas A numeric vector of lambda values where we will perform a spline fit on. Lambda values typically are in the range of 1e-12 to 1e-2 for interpolationg and for limiting lambda form respectively.
#' @param method Method for the projection index search. We always recommend using the default which is the grid approach for 1-dimeensional template manifolds. The alternative is "optimize" based on built-in R routines. See function project_optimize.
#' @param gridSize the number of data points to define the partition of the grid search when finding optimal projection indices. See section C.2 "computing projection indices" in the paper
#' @param min_iter the minimum number of total iterations of PME. Default is 1 for speed, but choosing a higher number is useful if you want to visualize the cost function.
#' @param max_iter the maximum number of iterations of PME. Prevents the algorithm from stalling.
#' @param err_tol error tolerance for convergence of the algorithm. Lower values result in longer runtimes, and more likely to hit the maximum iteration cap.
#' @param init_iso Whether or not to use the ISOMAP initialization. This should be set to false if the data satisfies a star-shape property (e.g. circles, ellipses, star-shapes, etc). If not, set to true (e.g., moon shapes, cashew shapes, etc.). See section C.1 Initialization in the paper for more details.
#' @param optimize_lambda whether or not to perform the method of section 5, Tuning Parameter Selection. If true, additional values will be returned, being the coefficient of variation at each lambda, along with the mean residual and standard deviation residual used to calculate it.
#' @param cl A a makeCluster object for parallelization. For most cases, should be left null so that a new cluster is made upon runtime.
#' @param lvp_grid The number denoting the size of the grid used for kernel smoothing for the tuning parameter selection. If left as NULL, it will choose a grid of 10 times the sample size.
#' @param parallel whether or not to run the algorithm parallel. The algorithm will be parallelized with respect to lambda and is recommended if your lambdas vector contains many values. Not recommended for only a handful of lamda values unless your dataset is large.
#' @param k_init The starting value for a search of the ISOMAP nearest neighbor parameter k. Our method is to find the smallest such k such that the ISOMAP graph is connected, as detailed in the iso_try function. The default is 1 to find the smallest possible.
#' @param approx whether or not to use the approximate spline method. We recommend to set this to be TRUE for very large lambda vectors and or very large sample sizes.
#' @export
#' @importFrom foreach %dopar% %do%
pme_1d_periodic <- function(data, lambdas, method = "grid", gridSize = 1000,  min_iter=1, max_iter = 10, err_tol=0.01, k_init = 1,init_iso = F, optimize_lambda=F,
                                  parallel=F,cl=NULL,lvp_grid=NULL, approx=F)
{
  n <- nrow(data)
  X <- data
  D <- ncol(data)
  
  
  lvp_grid <- ifelse(is.null(lvp_grid), 10*n, lvp_grid)
  
  # Initialize step: project onto unit circle. This is independent of the lambda choice.
  tcircle <- NULL
  if(D == 2 & !init_iso){
    tcircle <- project_circle(X)
  }else{
    iso2t <- iso_try(X, k_init, 2)
    
    iso2 <- iso2t$iso
    k <- iso2t$k
    print(paste0("k is ", toString(k)))
    
    tcircle <- project_circle(iso2$points)
  }
  
  
  # Initialize step: project onto unit circle. This is independent of the lambda choice.
  tinit <- NULL
  
  # Rule of thumb for choosing iso. Smallest possible without error
  isores <- iso_try(X, k_init)
  iso <- isores$iso
  k <- isores$k
  print(paste0("k is ", toString(k)))
  
  tinit <- normalize(iso$points)
  
  
  miniproj <- function(Xraw, newSpline){
    Xt<-NULL
    if(method == "grid"){
      tt <- project_grid(Xraw, newSpline, gridSize)
      #print(tt)
    } else if(method == "optimize"){
      tt<- project_optimize(Xraw,newSpline,periodic=T)
    }
    tt
  }
  
  spline1d_periodic_func <- ifelse(approx, spline1d_periodic_approx,
                                   spline1d_periodic)
  # Let fk denote kth spline fit, and pi_k the kth projection.
  # fk is based on pi-(k-1), and pi_k is based on fk.
  # Repeat whole scheme for each lambda
  
  # Setup parallel cluster
  # Setup cluster if not provided
  stop_after <- FALSE
  if (is.null(cl) && parallel) {
    numCores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(numCores)
    doParallel::registerDoParallel(cl)
    stop_after <- TRUE
  }
  # run parallel? T or F
  operator <- if(parallel) `%dopar%` else `%do%` #thanks chatgpt
  
  # Main algorithm loop
  results <-operator( foreach::foreach(i = 1:length(lambdas), 
                              .packages = c("stats", "pracma"),
                              .export = c("spline1d_periodic","project_grid","project_circle", "project_optimize",
                                          "spline_error","krt","krst", "local_var_1d","var_het")), {
                                            
                                            print("init start")
                                            lambda <- lambdas[i]
                                            #splines <- vector(mode="list",length=max_iter)
                                            alphas <-rep(0,max_iter)
                                            fidelities <- rep(0, max_iter)
                                            
                                            # Setup first spline fit
                                            spl <- spline1d_periodic_func(X,tcircle, lambdas[i])
                                            alphas[1] <- attr(spl, "alpha")
                                            
                                            
                                            # Setup first projection 
                                            # tt <- list() # list of t values
                                            t_proj<- miniproj(X,spl) %% 1 # evaluate first projection set
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
                                              spl <- spline1d_periodic_func(X, t_proj, lambdas[i])
                                              alphas[iter] <- attr(spl, "alpha")
                                              
                                              # Calculate new projection
                                              t_proj <- miniproj(X,spl) %% 1 # evaluate first projection set
                                              fidelities[iter] <- spline_error(X,t_proj, spl)
                                              
                                              # Get newest error
                                              #Diter[iter] <- spline_error(X,tt[[iter]], splines[[iter]])
                                              
                                              # Compare old error
                                              err <-abs((fidelities[iter] - fidelities[iter-1])/max(fidelities[iter-1], 1e-12))
                                            }
                                            
                                            resid_sq <- rowSums( (X - spl(t_proj))^2)
                                            
                                            # Compute variance heterogeneity
                                            vh <- NULL
                                            if(optimize_lambda){
                                              lvp <- local_var_1d(t_proj, resid_sq,grid_eval=lvp_grid, periodic=T)$sigma2_hat
                                              vh <- var_het(lvp) # get coefficient of variation
                                            }
                                            
                                            
                                            list(spline = spl, X = X, tproj = t_proj, lambda = lambdas[i],
                                                 resid_sq = resid_sq, cv = vh$cv,alphas=alphas[1:iter],
                                                 mean_sq = fidelities[1:iter],cost_functional = lambdas[i]*alphas[1:iter] + fidelities[1:iter],
                                                 vh_mean = vh$mean,
                                                 vh_sd = vh$sd)
                                          })
  # Aggregate results and find minimizing spline
  spline_list <- lapply(results, function(res) res$spline)
  lambda_list <- lambdas
  proj_list <- lapply(results, function(res) res$tproj)
  alpha_list <- lapply(results, function(res) res$alphas)
  mean_sq_list <- lapply(results, function(res) res$mean_sq)
  cost_functional_list <- lapply(results, function(res) res$cost_functional)
  
  list_return <- list(spline_list = spline_list, lambda_list = lambda_list, proj_list=proj_list,
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
    if (stop_after) parallel::stopCluster(cl)
  }, add = TRUE)
  
  
  list_return$template <- "S1"
  list_return$D <- D
  class(list_return) <- "cpme"
  
  return(list_return)
}
