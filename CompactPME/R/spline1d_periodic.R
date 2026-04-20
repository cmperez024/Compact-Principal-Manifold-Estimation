# Generated from create-CompactPME.Rmd: do not edit by hand

#' Fit spline
#' 
#' @param X data
#' @param t projection indices
#' @param lambda smoothing parameter
#' @export
spline1d_periodic <- function(X,t, lambda) {
  n <- nrow(X)
  #t <- as.vector(t)
  t<-as.vector(t)%% 1# normalize to 0 1. CHANGE MENG CHANGE MENG
  # ambient dimension is D
  
  # define matrices. m = 2 corresponds to 2nd deriv penaltty
  m <- 2
  r <- 0:m
  Tmat <- ones(n,1)
  #K <-  ((-1)^(m - 1)) *combn(1:n, n, Vectorize(function(j, k) krst(2 * m, t[j], t[k])))[,,1]
  K <- (-1)^(m-1) * outer(t, t, function(s, t) krst(2*m, s, t))
  M <- K + n * lambda * eye(n)
  
  # Theta factor (multiplies y) same for x and y coordinate. Calculate once.
  #Minv <- solve(M)
  # K is symmetric positive definite by propertie of rk. M is also symmetric positive definite
  # chol2inv chol M is 2x faster than normal solve
  Minv <- chol2inv(chol(M))
  
  # Precompute coefficients
  theta_coeffs <- ( solve(t(Tmat) %*% Minv %*% Tmat) %*% t(Tmat) %*% Minv) %*% X           # (p x D)
  alpha_coeffs <- Minv %*% (X - Tmat %*% theta_coeffs)  # (n x D)
  
  # calculate periodic function each dimension. Help of chatgpt ===
  fx <- function(x) {
    # Compute kernel matrix once
    Fmat <- outer(x, t, function(xi, tj) krst(2*m, xi, tj))  # length(x) x n
    # Combine all dimensions at once
    (-1)^(m-1) * Fmat %*% alpha_coeffs + matrix(rep(theta_coeffs, each=length(x)), nrow=length(x))
    # end chatgpt ===
  }
  
  attr(fx, "alpha") <- sum(diag(t(alpha_coeffs) %*% K %*% alpha_coeffs))
  
  return(fx)
  # return function list for each dim
}

#' Approximate spline for S1
#' @param X NxD data matrix
#' @param t projection indices
#' @param lambda smoothing parameter value
#' @param k number of basis functions for approximation
#' @param rescale_lam whether or not to rescale lambda so that the lambda choices for the approximate method correspond roughly to the exact method. Setting to T tries to make the scale more similar to the exact method, so identical lambda across the methods yield similar results
#' @export
spline1d_periodic_approx <- function(X,t,lambda, k=60,rescale_lam=F){
  
  # Get spline fits for arbitrary ambient dimension
  D <- ncol(X)
  N <- nrow(X)
  fits <- vector(mode="list", length=D)
  
  #t <- as.vector(t)
  t<-as.vector(t)%% 1 # normalize to 0 1. CHANGE MENG CHANGE MENG
  
  # scale lambda to our typical domain
  mylambda <- lambda* (10^9)^rescale_lam
  fits <- apply(X,2, function(data_col){
    
    mgcv::gam(data_col ~ s(t, bs = "cc", k=k),
              knots = list(t=c(0,1)),
              sp = mylambda,
              method = "REML",
              #norescale=T
              )
  })
  
  
  # Calculate penalty
  # Return penalty
  penalty <- 0
  for(i in seq_along(fits)){
    fiti <- fits[[i]]
    S <- fiti$smooth[[1]]
    # Get  penalized part of coefficient vector
    c <- fiti$coefficients[S$first.para:S$last.para]
    penalty <- penalty + t(c) %*% S$S[[1]] %*% c
  }
  
  
  spline_fun <- function(tnew) {
    eval <- lapply(fits, function(spline) predict(spline, newdata=data.frame(t=tnew),type="response"))
    do.call(cbind, eval)
  }
  
  attr(spline_fun, "alpha") <- penalty 
  spline_fun
  
}
