# Generated from create-CompactPME.Rmd: do not edit by hand

#' # The spline fit for the 1d Interval
#' 
#' @param X the NxD data matrix where N is the number of observations and D is either 2 or 3.
#' @param t the set of projection indices associated with the observations of X. Should be a vector of length N.
#' @param lambda the regularization parameter. Should be a real number greater than 0, typically in the magnitude of 1e-12 (interpolating) to 1e-2 (corresponds to lambda=infinity case in paper).
#' @return Returns the spline fit as a function and alpha value as an attribute.
#' @export
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

#' # The approximate spline fit for the 1d Interval
#' 
#' @param X the NxD data matrix where N is the number of observations and D is either 2 or 3.
#' @param t the set of projection indices associated with the observations of X. Should be a vector of length N.
#' @param lambda the regularization parameter.
#' @return The spline fit as a function and alpha value as an attribute.
#' @export
spline1d_interval_approx <- function(X,t,lambda){
  
  # Get spline fits for arbitrary ambient dimension
  D <- ncol(X)
  N <- nrow(X)
  fits <- vector(mode="list", length=D)
  
  fits <- apply(X,2, function(data_col){
    smooth.spline(t, data_col,lambda=lambda)
  })
  
  
  # Calculate penalty matrix for use in the algorithm to visualize cost functional
  # We utilize the fact that the approximate method is based on B-spline expansion
  # Return penalty
  penalty <- 0
  for(i in seq_along(fits)){
    fiti <- fits[[i]]
    c <- fiti$fit$coef
    knotvals <- fiti$fit$knot
    
    t_grid <- seq(min(knotvals), max(knotvals), length.out=10*N)  # fine grid
    delta <- diff(t_grid)[1]  # step size for numerical integration
    
    B_dd <- splines::splineDesign(
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
