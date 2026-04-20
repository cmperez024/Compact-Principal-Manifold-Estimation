# Generated from create-CompactPME.Rmd: do not edit by hand

#' An isomap method to find the smallest k so that the graph is not fragmented
#' @param X the NxD data matrix
#' @param k the starting value for the k parameter search
#' @param d The dimension we wish to get values for. Use 1 for the line, 2 for S1, 3 for S2
#' @return A list consisting of the `isomap` result and the optimal `k` value.
iso_try <- function(X, k=1, d=1){
  k<-k
  repeat{
    output <- try(vegan::isomap(dist(X),ndim=d,k=k),silent = T)
    if(!inherits(output, "try-error")) return(list(iso=output, k = k))
    
    k <- k+1
    
  }
}

#' Normalizing for the 0,1 interval
#' @param t vector of projection indices
#' @return The normalized values
normalize <- function(t){
  a <- min(t)
  b <- max(t)
  (t-a)/(b-a)
}

#' Computes residuals. Also works for the 2d case
#' @param Xdata NxD data matrix
#' @param t vector of projection indices
#' @param splineresult a spline fit as a function
#' @return The mean of squared residuals
spline_error <- function(Xdata, t, splineresult){
  mean((rowSums( (Xdata-splineresult(t))^2)))
}

