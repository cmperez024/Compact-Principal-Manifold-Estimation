# Generated from create-CompactPME.Rmd: do not edit by hand

#' Circular projection, for initialization puprposes
#' @param X A data `matrix` consisting of N rows and 2 columns
#' @return A vector of N values between 0 and 1 denoting the circular projection index
project_circle <- function(X){
  Xc <- scale(X, scale=F)
  (atan2(Xc[,2],Xc[,1])/(2*pi)) %% 1
}



#' Optimization over a grid. This is to find the best projection index on a 1d Manifold
#' @param Xraw A data `matrix`
#' @param res A spline fit function
#' @param gridSize The number of entries in the partition
#' @return A vector of projection indices onto the spline
project_grid <- function(Xraw, res, gridSize = 1000) {
  # Initialize grid
 # t_grid <- seq(0, 1, length.out = gridSize) UNCOMMENT THIS UNCOMMENT REMOVE THIS
  t_grid <- seq(0, 1, length.out = gridSize+1)[-(gridSize+1)]
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

#' Find projection index using an optimization routine. Grid method is preferred for global solutions
#' @param Xraw A data `matrix`
#' @param res A spline fit function
#' @param periodic Whether or not the spline is periodic. `T` for the circle, `F` for the interval
#' @return The projection indices
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


