# Generated from create-CompactPME.Rmd: do not edit by hand

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

