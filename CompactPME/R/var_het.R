# Generated from create-CompactPME.Rmd: do not edit by hand

#' Computes variance heterogneity given the results of local_var_1d. Also used for the 2d case.
#' @param sigma2_hat vector of smoothed MSE values
#' @return A list consisting of the mean, standard deviation, coefficient of variation, and spread.
var_het <- function(sigma2_hat) {
  mean_sigma2 <- mean(sigma2_hat)
  sd_sigma2   <- sd(sigma2_hat)
  cv          <- sd_sigma2 / (mean_sigma2 + 1e-12)  # coefficient of variation
  spread      <- max(sigma2_hat) - min(sigma2_hat)
  list(mean = mean_sigma2, sd = sd_sigma2, cv = cv, spread = spread)
}


#'  Compute smoothed residuals using kernel smoothing
#'  @param t projection indices
#'  @param resid_sq residual values
#'  @param grid_eval grid size for smoothing
#'  @param bandwidth bandwidth value. NULL for automatic selection
#'  @param kernel the kernel method for ksmooth
#'  @return A list consisting of the evaluation grid, the fitted residuals, and the bandwidth.
local_var_1d <- function(t, resid_sq, periodic=F,
                         grid_eval = 100,
                         bandwidth = NULL,
                         kernel = "normal") {
  ## periodic extension: t-1, t, t+1
  t_ext <- t %% 1 # CHANGE MENG added %% 1
  e2_ext <- resid_sq
  
  # If periodic (set periodic=T for circle template)
  if(periodic){
    t_ext  <- c(t_ext - 1, t_ext, t_ext + 1)
    e2_ext <-  rep(resid_sq, 3)
    print("Doing periodic")
  }
  
  ## ---------- automatic bandwidth selection ----------
  if (is.null(bandwidth)) {
    bw <- NA_real_
    
    if (requireNamespace("KernSmooth", quietly = TRUE)) {
      ord <- order(t)
      tx  <- t[ord]
      ey  <- resid_sq[ord]
      bw_try <- try(KernSmooth::dpill(x = tx, y = ey), silent = TRUE)
      
      if (!inherits(bw_try, "try-error") && is.finite(bw_try) && bw_try > 0) {
        bw <- bw_try
      }
    }
    
    if (!is.finite(bw) || bw <= 0) {
      bw <- stats::bw.SJ(t)
    }
    
    bandwidth <- bw
  }
  ## ---------------------------------------------------
 # grid <- seq(0, 1, length.out = grid_eval) CHANGE MENG
  grid <- seq(0, 1, length.out = grid_eval + 1)[-(grid_eval + 1)]  # [0,1)
  
  ks <- ksmooth(
    x        = t_ext,
    y        = e2_ext,
    kernel   = kernel,
    bandwidth = bandwidth,
    x.points = grid
  )
  sigma2_hat <- pmax(ks$y, 0)  # clamp to nonnegative
  list(s = grid, sigma2_hat = sigma2_hat, bandwidth = bandwidth)
}
