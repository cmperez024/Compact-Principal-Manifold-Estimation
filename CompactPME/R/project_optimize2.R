# Generated from create-CompactPME.Rmd: do not edit by hand

#' Computing projection indices for the spherical data case
#' @param X Nx3 data matrix, where N is the number of observations
#' @param f The function to optimize over. Usually a spline fit
#' @param init An N-vector consisting of initial values for the optimizing search
 project_optimize2 <- function(X, f, init) {
    param_return <- matrix(0, nrow = nrow(X), ncol = 2)
    
    init_sph <- pracma::cart2sph(init)[,-3]
    
    objective <- function(param,Xi) {
      input <- matrix(param, ncol=2)
      input_cart <- pracma::sph2cart(cbind(input, 1))
      sum((Xi - f(input_cart)  )^2  )
    }
    
    for (i in 1:nrow(X)) {
      nl <- nlm(objective, init_sph[i,], Xi=X[i,])$estimate
      param_return[i,] <- nl
    }
    pracma::sph2cart(cbind(param_return, 1))
 }


#' Sampling points on fibonacci sphere. Adapted from CoolTools library. Taken from github to avoid loading unnecessary functions
#' @param n sample size
#' @param r radius
fibonaccisphere = function(n=1000, r=1) {
  
  if (n<1 | round(n)!=n) stop('n must be a positive integer')
  
  goldenratio = (1+sqrt(5))/2
  i = seq(n)-0.5
  z = 1-2*i/n # z-coordinate for unit sphere
  theta = acos(pmax(-1,pmin(1,z))) # polar angle
  phi = (2*pi*i/goldenratio)%%(2*pi)
  
  x = r*sin(theta)*cos(phi)
  y = r*sin(theta)*sin(phi)
  z = r*z
  out = cbind(x=x,y=y,z=z)
  
  return(out)
  
}
