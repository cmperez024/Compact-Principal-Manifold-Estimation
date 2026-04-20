# Generated from create-CompactPME.Rmd: do not edit by hand

#'  qm helper function
#'  @param z The scalar value of some inner product
qm <- function(z){
    # if z>=1, issue
    # gpt suggestion. z forced to be between -1 and 1e-12
    z <- pmin(pmax(z,-1), 1-1e-12)
    zinv <- 1-z
    0.5*(log(1+sqrt(2/zinv))*( 3*zinv^2 - 2*zinv)-
           12*(zinv/2)^(1.5) + 6*(zinv/2)+1)
}

#' Another helper function for the sphere kernel
#' @param z The scalar value of some inner product
expr <- function(z) { 1/(4*pi) *(qm(z) - 1/3)}


#' Reproducing kernel for the sphere
#' @param v1 Vector
#' @param v2 Another vector
Rpp <- function(v1,v2){
  expr(sum(v1*v2))
} 
