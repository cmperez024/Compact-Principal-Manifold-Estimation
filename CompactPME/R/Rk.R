# Generated from create-CompactPME.Rmd: do not edit by hand

#' # The reproducing kernel for the interval
#' 
#' @param s eval location 1. Should be a vector of equal length to t. 
#' @param t eval location 2. Should be a vector of equal length to s.
#' @return The RK evaluation
Rk <- function(s,t){
  1/3*pmin(s,t)^3 - 1/2 *pmin(s,t)^2 * (s+t) + pmin(s,t)*s*t
}
