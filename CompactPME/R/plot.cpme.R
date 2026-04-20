# Generated from create-CompactPME.Rmd: do not edit by hand

#' Plotting function
#' @param fit A CPME object
#' @export
plot.cpme <- function(fit){
  
  template <- fit$template
  
  if(template == "interval" || template == "S1"){
    S <- seq(0,1,length.out=1000)
    plot_1d(fit, S)
  }
  else if(template =="S2"){
    S <- fibonaccisphere(1000)
    plot_2d(fit,S)
  }
  else{
    stop("Invalid template.")
  }
}

#' Plotting helper that calls the right function for 1d case
#' @param fit A CPME object
#' @param S The parameterization. A vector of numbers between 0 and 1 for interval and S1. kx3 vectors in R3 where k is the desired fidelity
#' @return A plotly object
plot_1d <- function(fit,S){
  
  data <- fit$dataset
  D <- fit$D
  
  # Start base plot with data points
  p <- plotly::yplot_ly()
  
  
  if(D == 2){
    p <- p |> plotly::add_trace(
      x = data[, 1], y = data[, 2],
      type = "scatter",
      mode = "markers",
      marker = list(size = 3, opacity = 0.7, color="black"),
      name = "Data points"
    )
  }else if(D==3){
    
    p <- p |> plotly::add_trace(
      x = data[, 1], y = data[, 2], z = data[, 3],
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 3, opacity = 0.7,color="black"),
      name = "Data points"
    )
    
  }
  
  # Loop through each fit and add its corresponding line
  for (i in seq_along(fit$spline_list)) {
    spl <- fit$spline_list[[i]]
    splinevis <- spl(S)
    lambda <- fit$lambda_list[i]
    label <- paste0("λ = ", formatC( fit$lambda_list[i],format="e",digits=2))
    
    if(D == 2){
      p <- p |> plotly::add_trace(
        x = splinevis[, 1], y = splinevis[, 2],
        type = "scatter",
        mode = "lines",
        marker = list(size = 2, opacity = 0.7),
        name = label
      )
    }else if(D==3){
      
      p <- p |> plotly::add_trace(
        x = splinevis[, 1], y = splinevis[, 2], z = splinevis[, 3],
        type = "scatter3d",
        mode = "lines",
        marker = list(size = 2, opacity = 0.7),
        name = label
      )
      
    }
  }
  
  return(p)
  
}


#' Helper for generating a surface
#' @param f A spline function
#' @return A `list` of evaluation matrices
genSurface <- function(f){
  M <- 50
  thetas <- seq(from = 0, to = 2*pi, length.out=M)
  phis <- seq(from =-pi/2, to = pi/2, length.out=M)
  rr <- 1
  grid <- expand.grid(thetas,phis,1)
  cartgrid <- t(apply(grid, 1, pracma::sph2cart))
  
  # Evaluate manifold
  fvalues <- f(cartgrid)
  
  Xmat <- matrix(fvalues[,1], nrow = M, ncol = M, byrow = TRUE)
  Ymat <- matrix(fvalues[,2], nrow = M, ncol = M, byrow = TRUE)
  Zmat <- matrix(fvalues[,3], nrow = M, ncol = M, byrow = TRUE)
  
  list(X=Xmat,Y=Ymat,Z=Zmat)
}


#' Helper for plotting a 2d surface
#' @param fit A CPME object
#' @param S a kx3 `matrix` consisting of k-many vectors in R3, denoting the evaluation points of the spline
#' @return Returns a `list` containing plotly objects
plot_2d <- function(fit,S){
  
  data <- fit$dataset
  plot_list <- vector("list", length=length(fit$spline_list))
 
  for (i in seq_along(fit$spline_list)) {
    
    # Get spline
    spl <- fit$spline_list[[i]]
    
    # Get surface
    surf <- genSurface(spl)
    
    lambda <- fit$lambda_list[i]
    label <- paste0("λ = ", formatC( fit$lambda_list[i],format="e",digits=2))
    
     plot_list[[i]] <- plotly::plot_ly() %>%
    # Scatter points trace
    plotly::add_markers(x = data[, 1],y = data[, 2], z = data[, 3],
              type = "scatter3d",
              mode = "markers",
              marker = list(size = 2,opacity=0.7, color ="black"),
              name = "Data points"
    ) %>%
    # Overlay line trace
    plotly::add_surface(x=surf[[1]],y=surf[[2]], z=surf[[3]], opacity = 0.3,showscale=F,colorscale="Viridis")|>
    plotly::layout(title = paste0("λ = ", formatC( fit$lambda_list[i],format="e",digits=2)),
           scene = list(
             xaxis = list(title = "X"),
             yaxis = list(title = "Y"),
             zaxis = list(title = "Z")
           ))
  }

  plot_list
  
}
