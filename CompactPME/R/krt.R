# Generated from create-CompactPME.Rmd: do not edit by hand

#' Helper function for the kernel
#' 
#' @param r degree of the kernel
#' @param t where we evaluate
krt <- function(r, t) {
  tfrac <- t - floor(t)
  bernoulli(r, tfrac) / factorial(r)
}

#' The Reproducing kernel for S1
#' 
#' @param r degree
#' @param s eval location 1
#' @param t eval location 2
krst <- function(r, s, t) {
  krt(r, s - t)
}
