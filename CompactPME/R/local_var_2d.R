# Generated from create-CompactPME.Rmd: do not edit by hand

#'  Compute smoothed residuals using kernel smoothing, 2d CASE
#'  @param Z projection indices (vectors in R3)
#'  @param resid_sq residual values
#'  @param eval_size grid size for smoothing
#'  @param bw_grid interval size for searching for optimal bandwidth
#'  @return A list consisting of the evaluation grid, the fitted residuals, and the bandwidth.
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


#' Circular distances between a single vector and every other vector in a matrix
#' @param x_new A vector denoting a new data point
#' @param x_data The dataset
#' @returns A vector of distances between the vector to each row of the matrix using a circular geodesic
d_circ_vec <- function(x_new, x_data) {
  # x_new: vector (length d)
  # x_data: matrix n_train x d
  dots <- as.vector(x_data %*% x_new)  # fast matrix multiply in R
  pmax(pmin(dots,1),-1) |> acos()
}

#' Implements the pelletier kernel using Gaussian kernel
#' @param x_news New x values to evaluate at
#' @param x_data The x data (projection indices)
#' @param y_data The y data (residuals)
#' @returns Estimated y values at the `x_news` values.
#' @references
#' Pelletier, B. (2006). Non-parametric regression estimation on closed Riemannian manifolds. Journal of Nonparametric Statistics, 18(1), 57-67.
pelletier_kernel <- function(x_news, x_data, y_data, h) {
  sapply(1:nrow(x_news), function(i){
    rho <- d_circ_vec(x_news[i,], x_data)
    k <- dnorm(rho / h)*ifelse(rho<1e-16, 1, rho/sin(rho))
    sum(k * y_data) / sum(k)
  })
}

