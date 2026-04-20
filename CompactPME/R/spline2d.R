# Generated from create-CompactPME.Rmd: do not edit by hand

#' Spline on the sphere
#' @param data An Nx3 data matrix
#' @param Tt An Nx3 projection index matrix. Each observation (row) is a vector parameterizing the location on the sphere S2 as opposed to longitude and latitude.
#' @param lambda The smoothing parameter. Real number greater than 0
#' @export
   spline2d <- function(data,Tt, lambda){
    n <- nrow(data)
  #Form R kernel matrix (chatgpt)
    # Normalize rows of data
    #Znorm <- data / sqrt(rowSums(data^2))  # n x 3
    #Tt <- scale(Tt,scale=F)
    Znorm <- Tt / sqrt(rowSums(Tt^2))
    # Compute all pairwise cos(gamma) = dot products
    cosGamma <- Znorm %*% t(Znorm)  # n x n
    
    #View(cosGamma)
    # Apply kernel function to all pairs at once
    Rmat2 <- matrix(0, nrow=n, ncol=n)
    # Diagonal: limit value
    diag(Rmat2) <- 1/(24*pi)
    
    # Off-diagonal: apply expr to all non-diagonal entries
    diag(cosGamma) <- NA  # ignore diagonal
    Rmat2[!is.na(cosGamma)] <- expr(cosGamma[!is.na(cosGamma)])
    
    #============== Precompute matrices
    m <- 2
    r <- 0:m
    Tmat <- ones(n,1)
    M <- Rmat2 + n * lambda * eye(n)
    Minv <- solve(M)
    
    # constant part
    const_factor <- solve( t(Tmat) %*% Minv %*% Tmat) %*% t(Tmat) %*% Minv %*% data
    coeffs <- Minv %*% (data - Tmat %*% const_factor) 
    
    
    #================== Return spline function
    f_eval_fast <- function(Pnew) {
      Pnew <- matrix(Pnew,ncol=3)
      
      # normalize rows
      Pnew_norm <- Pnew / sqrt(rowSums(Pnew^2))
      Z_norm <- Tt / sqrt(rowSums(Tt^2))
      
      # compute all pairwise dot products (m x n)
      cosGamma <- Pnew_norm %*% t(Z_norm)
      
      # apply kernel function elementwise
      # cosGamms is all dot products. We evaluate R kernel separately whether
      # cosGamma is 1 or not
      Rnew <- cosGamma
      oneloc <- which(cosGamma>=1, arr.ind=T)
      cloc <- which(cosGamma < 1, arr.ind=T)
      Rnew[oneloc] <- 1/(24*pi)
      Rnew[cloc] <- expr(Rnew[cloc])
      #new
      
      #Rnew <- expr(cosGamma)
      # evaluate spline
      Rnew %*% coeffs + matrix(rep(const_factor, each = nrow(Pnew)), nrow = nrow(Pnew))
    }
    
    # add penalty attribute
    attr(f_eval_fast, "alpha") <-sum(diag(t(coeffs) %*% Rmat2 %*% coeffs))
    # Return function
    f_eval_fast
    
    #return(Rmat2)
  }
