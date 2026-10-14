# Resolve paths from this file, including when invoked outside the project.
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_file)) stop("Run this script with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_file[1L]),
                                     winslash = "/", mustWork = TRUE))
figure_dir <- Sys.getenv("PME_OUTPUT_DIR", unset = file.path(script_dir, ".."))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- normalizePath(figure_dir, mustWork = TRUE)
repo_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."),
                           mustWork = TRUE)
local_library <- file.path(repo_root, "packages")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
setwd(figure_dir)
set.seed(1)
for (folder in c("results_1d", "results_sph", "job_figs_1d",
                  "panel_flower", "panel_cashew", "sim_results")) {
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
}

# Reproduce the appendix large-lambda convergence figure only.
# All three studies are computed here; no other simulation files are read.
required <- c("CompactPME", "vegan", "fdasrvf", "ggplot2", "patchwork")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " first.")
}
library(ggplot2)
library(patchwork)
isomap <- vegan::isomap
spline1d_periodic <- CompactPME::spline1d_periodic
spline1d_interval <- CompactPME::spline1d_interval
spline2d <- CompactPME::spline2d
project_circle <- CompactPME:::project_circle
project_grid <- CompactPME:::project_grid
project_optimize2 <- CompactPME:::project_optimize2
spline_error <- CompactPME:::spline_error
iso_try <- CompactPME:::iso_try
normalize <- CompactPME:::normalize
fibonaccisphere <- CompactPME:::fibonaccisphere
cart2sph <- pracma::cart2sph
sph2cart <- pracma::sph2cart
L2distmean <- function(sample_f, sample_mean) {
  sqrt(sum(sweep(sample_f, 2, sample_mean, "-")^2))
}
geodesicDistance <- function(sample, true, scale = FALSE) {
  fdasrvf::calc_shape_dist(t(sample), t(true), mode = "O",
                          scale = scale, rotation = FALSE)$d
}

rotation2 = function(angle) {
  c = cos(angle)
  s = sin(angle)
  return(rbind(c(c,-s),c(s,c)))
  
}

sine1d2D <- function(n, rotation=0, period=2*pi, noise=0, tval){
  t0 <- if(missing(tval)) sort(runif(n)) else tval
  
  list(X=cbind(t0, sin(period*t0)) %*% t(rotation2(rotation)) + rnorm(2*n, sd=noise),
       t = t0)
}

flower1d2D <- function(n, petals, noise, t){
  t0 <- if(missing(t)) sort(runif(n)) else t
  
  r <- (1 + 0.3 * sin(petals *2*pi* t0))
  return(list(t=t0, X = cbind(r * cos(2*pi*t0),  r * sin(2*pi*t0))+rnorm(2*n, sd=noise)))
}

flower2d3D <- function(N = 300,
                       r0 = 1, a = 0.3, petals = 6, b = 0.5, noise = 0.01, t) {
  if(!missing(t)){
    coords <- cart2sph(t)
    theta <- coords[,1]
    phi <- coords[,2]
    z <- t[,3]
  }else{
    theta <- runif(N, 0, 2*pi)
    z <- runif(N, -1, 1)
    phi <- acos(z)
  }
  
  r_xy <- r0 * (1 + a * cos(petals * theta))
  
  x <- r_xy * sqrt(1 - z^2) * cos(theta)
  y <- r_xy * sqrt(1 - z^2) * sin(theta)
  z <- b * z
  
  pts <- cbind(x, y, z) + rnorm(3*N, sd = noise)
  return(list(t = sph2cart(cbind(theta,phi,1)), X=pts))
}

pmeClosed_1d_periodic_simulation <- function(shape_type,D, Ns,noise, lambdas, gridSize = 1000, max_iter = 20, err_tol=0.01, iso_k = 10,init_iso = F, L2dist_size = 1000, parallel=F,cl=NULL){
  X <- vector(mode="list", length=length(Ns))
  trues <- vector(mode="list", length=length(Ns))
  tcircles <- vector(mode="list", length=length(Ns))
  
  
  # the grid to evaluate at. this is also where we draw the data from.
  L2dist_grid <- seq(from=0,to=1-1e-6, length.out=L2dist_size)
  
  
  # Initialize data and true sample data
  if(shape_type == "ellipse" & D==2){
    
    for(j in seq_along(Ns)){
      
      X[[j]] <- circle1d2D(Ns[j], 1, 2, c(0,0), noise)$X
      trues[[j]] <- colMeans(X[[j]])#circle1d2D(L2dist_size, 1, 2, c(0,0), 0, t=L2dist_grid)$X
    }
    
  }else if(shape_type =="flower" & D==2){
    
    for(j in seq_along(Ns)){
      X[[j]] <- flower1d2D(Ns[j], 5, noise)$X
      trues[[j]] <- colMeans(X[[j]])##flower1d2D(L2dist_size,5, 0, t=L2dist_grid)$X
    }  
    
  }
  
  # Get initial projection indices for each dataset
  for(i in 1:length(Ns)){
    if(D == 2 & !init_iso){
      tcircles[[i]] <- project_circle(X[[i]])
    }else{
      iso2 <- isomap(dist(X[[i]]), ndim=2, k=iso_k)
      tcircles[[i]] <- project_circle(iso2$points)
    }
      
  }
  
  
  miniproj <- function(Xraw, newSpline){
    project_grid(Xraw, newSpline, gridSize)
  }
  
  # Setup parallel cluster
  stop_after <- FALSE
  if (is.null(cl) && parallel) {
    numCores <- parallel::detectCores() - 1
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    stop_after <- TRUE
  }
  
  # run parallel? T or F
  # Lambdas are evaluated sequentially with lapply.
  
  
  results <- vector(mode="list", length = length(Ns))
  
  # edit name to index by N
  names(results) <- paste0("N_",Ns)
  
  for(j in seq_along(Ns)){
    
    
    # Main algorithm loop
    results[[j]] <-lapply(seq_along(lambdas), function(i) {
                                                   lambda <- lambdas[i]
                                                   splines <- list()
                                                   L2distances <- rep(NA, max_iter+1)
                                                   geodesicdistance <- rep(NA, max_iter+1)
                                                   L2index <- 1
                                                   
                                                   # Setup first spline fit
                                                   splines[[1]] <- spline1d_periodic(X[[j]],tcircles[[j]], lambdas[i])
                                                   
                                                   # Get first true manifold error
                                                   L2distances[1] <- L2distmean(splines[[1]](L2dist_grid), trues[[j]])
                                                   #geodesicdistance[1] <- geodesicDistance(splines[[1]](L2dist_grid), trues[[j]])
                                                   L2index <- L2index + 1
                                                   
                                                   # Setup first projection 
                                                   tt <- list() # list of t values
                                                   tt[[1]] <- miniproj(X[[j]],splines[[1]]) # evaluate first projection set
                                                   
                                                   # Repeat until error small
                                                   err <- err_tol + 1
                                                   iter <- 1
                                                   Diter <- rep(Inf, max_iter)
                                                   # Initial error
                                                   Diter[iter] <- spline_error(X[[j]], tt[[iter]], splines[[iter]])
                                                   
                                                   # Iteriate until error convergences
                                                   while(err > err_tol &  iter < max_iter){
                                                     iter <- iter+1
                                                     
                                                     #  Get newest spline fit
                                                     splines[[iter]] <- spline1d_periodic(X[[j]], tt[[iter-1]], lambdas[i])
                                                     # Get  true manifold error
                                                     L2distances[L2index] <- L2distmean(splines[[iter]](L2dist_grid), trues[[j]])
                                                     #geodesicdistance[L2index] <- geodesicDistance(splines[[iter]](L2dist_grid), trues[[j]])
                                                     L2index <- L2index + 1
                                                     
                                                     # Calculate new projection
                                                     tt[[iter]] <- miniproj(X[[j]],splines[[iter]]) # evaluate first projection set
                                                     
                                                     # Get newest error
                                                     Diter[iter] <- spline_error(X[[j]],tt[[iter]], splines[[iter]])
                                                     
                                                     # Compare old error
                                                     err <-abs((Diter[iter] - Diter[iter-1])/Diter[iter-1])
                                                   }
                                               
                                                   lambda_res <- list(
                                                     lambda = lambdas[i],
                                                    # spline_sequence = list(
                                                       splines = splines,
                                                       projections = tt,
                                                       max_iter = iter,
                                                       L2 = L2distances[1:iter]
                                                       #geodesic = geodesicdistance[1:iter]
                                                    # ),
                                                     #diagnostics = list(
                                                     #  L2 = L2distances[1:iter],
                                                    #   geodesic = geodesicdistance[1:iter]
                                                   #  )
                                                   )
                                                
                                                 })
    
    # Rename results j
    names(results[[j]]) <- paste0("lambda_", lambdas)
    results[[j]]$data <- X[[j]]
  }
  
    on.exit({
    if (stop_after) stopCluster(cl)
  }, add = TRUE)

  results
}

pmeClosed_1d_interval_simulation <- function(shape_type, D, Ns, noise, lambdas, gridSize = 1000, max_iter = 10, err_tol=0.01, cl=NULL,parallel=F, k_init = 1, L2dist_size = 1000)
{
  # Begin actual code 
  X <- vector(mode="list", length=length(Ns))
  trues <- vector(mode="list", length=length(Ns))
  
  tinits <- vector(mode="list", length=length(Ns))
  
  # L2 dist initialize
  L2dist <- function(sample_f, true_f){
    sqrt(sum((sample_f-true_f)^2))
  }
  
  # the grid to evaluate at. this is also where we draw the data from.
  L2dist_grid <- seq(from=0,to=1-1e-6, length.out=L2dist_size)
  

  if(shape_type == "sine" & D==2){
    
    for(j in seq_along(Ns)){
      
      X[[j]] <- sine1d2D(Ns[j], rotation=pi/6,period=2*pi, noise=noise)$X
      trues[[j]] <- sine1d2D(L2dist_size, rotation=pi/6,period=2*pi, noise=noise)$X
    }
  }else if(shape_type =="PCA" & D==2){
    
    for(j in seq_along(Ns)){
      X[[j]] <- sine1d2D(Ns[j], rotation=pi/6,period=2*pi, noise=noise)$X
      
      # Get PCA line
      pca <- prcomp(X[[j]], center=T)
      pca1 <- pca$rotation[,1]
      slope <- pca1[2]/pca1[1]
      mu <- colMeans(X[[j]])
      intr <- mu[2] - slope*mu[1]
      
      trues[[j]] <- cbind(L2dist_grid, slope*L2dist_grid + intr)
    }  
    
  }
  use_scale <- shape_type=="PCA"
  
  # Initialize first projections 
  k_init<-1
  
  for(i in 1:length(Ns)){
    isores <- iso_try(X[[i]], k_init)
    iso <- isores$iso
    k <- isores$k
    print(paste0("k is ", toString(k)))
    tinits[[i]] <-normalize(iso$points)
  }


  # Setup parallel cluster
  # Setup cluster if not provided
  stop_after <- FALSE
  if (is.null(cl) && parallel) {
    numCores <- parallel::detectCores() - 1
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    stop_after <- TRUE
  }
  
  
  miniproj <- function(Xraw, newSpline){
    project_grid(Xraw, newSpline, gridSize)
  }
  # run parallel? T or F
  # Lambdas are evaluated sequentially with lapply.
  
  
  
  results <- vector(mode="list", length = length(Ns))
  names(results) <- paste0("N_",Ns)
  
  for(j in seq_along(Ns)){
    
    
    
    # Main algorithm loop
    results[[j]] <-lapply(seq_along(lambdas), function(i) {
                                                   lambda <- lambdas[i]
                                                   splines <- list()
                                                   L2distances <- rep(NA, max_iter+1)
                                                   geodesicdistance <- rep(NA, max_iter+1)
                                                   L2index <- 1
                                                   
                                                   # Setup first spline fit
                                                   splines[[1]] <- spline1d_interval(X[[j]],tinits[[j]], lambdas[i])
                                                   
                                                   # Get first true manifold error
                                                   L2distances[1] <- L2dist(splines[[1]](L2dist_grid), trues[[j]])
                                                   geodesicdistance[1] <- geodesicDistance(splines[[1]](L2dist_grid), trues[[j]],use_scale)
                                                   L2index <- L2index + 1
                                                   
                                                   
                                                   
                                                   # Setup first projection 
                                                   tt <- list() # list of t values
                                                   tt[[1]] <-  normalize(miniproj(X[[j]],splines[[1]])) # evaluate first projection set
                                                   
                                                   # Repeat until error small
                                                   err <- err_tol + 1
                                                   iter <- 1
                                                   Diter <- rep(Inf, max_iter)
                                                   # Initial error
                                                   Diter[iter] <- spline_error(X[[j]], tt[[iter]], splines[[iter]])
                                                   
                                                   # Iteriate until error convergences
                                                   while(err > err_tol &  iter < max_iter){
                                                     iter <- iter+1
                                                     
                                                     #  Get newest spline fit
                                                     splines[[iter]] <- spline1d_interval(X[[j]], tt[[iter-1]], lambdas[i])
                                                   
                                                     # Get  true manifold error
                                                     L2distances[L2index] <- L2dist(splines[[iter]](L2dist_grid), trues[[j]])
                                                     geodesicdistance[L2index] <- geodesicDistance(splines[[iter]](L2dist_grid), trues[[j]],use_scale)
                                                     L2index <- L2index + 1
                                                     
                                                     # Calculate new projection
                                                     tt[[iter]] <-  normalize(miniproj(X[[j]],splines[[iter]]))
                                                     
                                                     # Get newest error
                                                     Diter[iter] <- spline_error(X[[j]],tt[[iter]], splines[[iter]])
                                                     
                                                     # Compare old error
                                                     err <-abs((Diter[iter] - Diter[iter-1])/Diter[iter-1])
                                                   }
                                                   
                                                   
                                                   lambda_res <- list(
                                                     lambda = lambdas[i],
                                                     # spline_sequence = list(
                                                     splines = splines,
                                                     projections = tt,
                                                     max_iter = iter,
                                                     L2 = L2distances[1:iter],
                                                     geodesic = geodesicdistance[1:iter]
                                                     # ),
                                                     #diagnostics = list(
                                                     #  L2 = L2distances[1:iter],
                                                     #   geodesic = geodesicdistance[1:iter]
                                                     #  )
                                                   )
                                                   
                                                 })
    
    # Rename results j
    names(results[[j]]) <- paste0("lambda_", lambdas)
    results[[j]]$data <- X[[j]]
  }
  
    on.exit({
    if (stop_after) stopCluster(cl)
  }, add = TRUE)

  results
}

pmeClosed_2d_simulation <- function(shape_type, Ns, noise, lambdas, gridSize = 1000, max_iter = 10, err_tol=0.01, cl=NULL,parallel=F, iso_k=10, L2dist_size = 1000){
  
  # Begin actual code 
  X <- vector(mode="list", length=length(Ns))
  trues <- vector(mode="list", length=length(Ns))
  
  tspheres <- vector(mode="list", length=length(Ns))
  
  # L2 dist initialize
  L2dist <- function(sample_f, true_f){
    sqrt(sum((sample_f-true_f)^2))
  }
  
  # the grid to evaluate at. this is also where we draw the data from.
  L2dist_grid <- fibonaccisphere(n=L2dist_size)
  
  
  if(shape_type == "ellipse"){
    
    for(j in seq_along(Ns)){
      
      X[[j]] <- ellipse2d3D(Ns[j], c(1,2,1), noise=noise)
      trues[[j]] <- colMeans(X[[j]])
    }
  }else if(shape_type =="flower"){
    
    for(j in seq_along(Ns)){
      
      X[[j]] <- flower2d3D(Ns[j], petals=5, noise=noise)$X
      trues[[j]] <-colMeans(X[[j]])
      
    } 
  }
  
  
  # Initialize first projections 
  
  for(i in 1:length(Ns)){
    
    if(shape_type != "moon"){
      datac <- scale(X[[i]],scale=F)
      tspheres[[i]] <- datac/sqrt(rowSums(datac^2))
    }
    else{
      iso3 <- isomap(dist(X[[i]]), ndim=3, k=iso_k)
      t3 <- scale(iso3$points,scale=F)
      spheres[[i]] <-t3 / sqrt(rowSums(t3^2))
    }
  }


  # Setup parallel cluster
  # Setup cluster if not provided
  stop_after <- FALSE
  if (is.null(cl) && parallel) {
    numCores <- parallel::detectCores() - 1
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    stop_after <- TRUE
  }
  
  
  miniproj <- function(Xraw, newSpline){
    project_grid(Xraw, newSpline, gridSize)
  }
  # run parallel? T or F
  # Lambdas are evaluated sequentially with lapply.
  
  
  
  results <- vector(mode="list", length = length(Ns))
  names(results) <- paste0("N_",Ns)
  
  for(j in seq_along(Ns)){
    
    
    
    # Main algorithm loop
    results[[j]] <-lapply(seq_along(lambdas), function(i) {
                                                   lambda <- lambdas[i]
                                                   splines <- list()
                                                   L2distances <- rep(NA, max_iter+1)
                                                   geodesicdistance <- rep(NA, max_iter+1)
                                                   L2index <- 1
                                                   
                                                   # Setup first spline fit
                                                   splines[[1]] <- spline2d(X[[j]],tspheres[[j]], lambdas[i])
                                                   
                                                   # Get first true manifold error
                                                   L2distances[1] <- L2distmean(splines[[1]](L2dist_grid), trues[[j]])
                                                   #geodesicdistance[1] <- geodesicDistance(splines[[1]](L2dist_grid), trues[[j]],use_scale)
                                                   L2index <- L2index + 1
                                                   
                                                  # Setup first projection 
                                                   tt <- list() # list of t values
                                                   tt[[1]] <-  project_optimize2(X[[j]], splines[[1]], tspheres[[j]]) # evaluate first projection set
                                                   
                                                   # Repeat until error small
                                                   err <- err_tol + 1
                                                   iter <- 1
                                                   Diter <- rep(Inf, max_iter)
                                                   # Initial error
                                                   Diter[iter] <- spline_error(X[[j]],tt[[iter]], splines[[iter]])
                                                   
                                                   # Iteriate until error convergences
                                                   while(err > err_tol &  iter < max_iter){
                                                     iter <- iter+1
                                                     
                                                     #  Get newest spline fit
                                                     splines[[iter]] <- spline2d(X[[j]], tt[[iter-1]], lambdas[i])
                                                   
                                                     # Get  true manifold error
                                                     L2distances[L2index] <- L2distmean(splines[[iter]](L2dist_grid), trues[[j]])
                                                     #geodesicdistance[L2index] <- geodesicDistance(splines[[iter]](L2dist_grid), trues[[j]],use_scale)
                                                     L2index <- L2index + 1
                                                     
                                                     # Calculate new projectio
                                                     tt[[iter]] <-  project_optimize2(X[[j]], splines[[iter]], tt[[iter-1]])
                                      
                                                     # Get newest error
                                                     Diter[iter] <- spline_error(X[[j]],tt[[iter]], splines[[iter]])
                                                     
                                                     # Compare old error
                                                     err <-abs((Diter[iter] - Diter[iter-1])/Diter[iter-1])
                                                   }
                                                   
                                                   
                                                   lambda_res <- list(
                                                     lambda = lambdas[i],
                                                     # spline_sequence = list(
                                                     splines = splines,
                                                     projections = tt,
                                                     max_iter = iter,
                                                     L2 = L2distances[1:iter]
                                                    # geodesic = geodesicdistance[1:iter]
                                                     # ),
                                                     #diagnostics = list(
                                                     #  L2 = L2distances[1:iter],
                                                     #   geodesic = geodesicdistance[1:iter]
                                                     #  )
                                                   )
                                                   
                                                 })
    
    # Rename results j
    
    names(results[[j]]) <- paste0("lambda_", lambdas)
    results[[j]]$data <- X[[j]]
  }
  results
}

smoke <- identical(Sys.getenv("PME_SMOKE_TEST"), "1")
Ns <- if (smoke) c(30L, 40L) else c(100L, 200L, 500L, 1000L)
noise <- 0.03
lambdas1 <- if (smoke) c(1e-6, 1) else 10^(-10:2)
lambdas3 <- if (smoke) c(1e-6, 1) else 10^(-10:0)
iterations <- if (smoke) 2L else 10L
set.seed(1)
results1 <- pmeClosed_1d_periodic_simulation(
  "flower", 2, Ns, noise, lambdas1, max_iter = iterations,
  L2dist_size = if (smoke) 40L else 300L, parallel = FALSE)
set.seed(1)
results2 <- pmeClosed_1d_interval_simulation(
  "PCA", 2, Ns, noise, lambdas1, max_iter = iterations,
  L2dist_size = if (smoke) 40L else 2000L, parallel = FALSE)
set.seed(1)
results3 <- pmeClosed_2d_simulation(
  "flower", Ns, noise, lambdas3, max_iter = iterations,
  L2dist_size = if (smoke) 40L else 1000L, parallel = FALSE)

to_frame <- function(results, lambdas, metric = "L2") {
  do.call(rbind, lapply(seq_along(Ns), function(j) {
    do.call(rbind, lapply(seq_along(lambdas), function(i) {
      data.frame(N = factor(Ns[j], levels = Ns), lambda = lambdas[i],
                 distance = tail(results[[j]][[i]][[metric]], 1L))
    }))
  }))
}
draw_panel <- function(data, title, label) {
  ggplot(data, aes(x = lambda, y = distance, col = N, linetype = N)) +
    geom_point(aes(shape = N), size = 2) +
    geom_line(linewidth = 0.8) + scale_x_log10() + theme_bw() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    xlab("λ") + ylab(label) + ggtitle(title) +
    theme(panel.grid = element_blank(), aspect.ratio = 1,
          legend.background = element_rect(fill = "white", color = "black"))
}
PCA1d <- draw_panel(to_frame(results2, lambdas1, "geodesic"),
                    "(i) Interval", "Geodesic Distance to PCA Line")
Sample1d <- draw_panel(to_frame(results1, lambdas1),
                       "(ii) Circle", "L2 Distance to Sample Mean")
Sample2d <- draw_panel(to_frame(results3, lambdas3),
                       "(iii) 2-Sphere", "L2 Distance to Sample Mean")
myplt <- (PCA1d + Sample1d + Sample2d + plot_layout(guides = "collect")) &
  theme(legend.position = "bottom")
ggsave(file.path(figure_dir, "dist_lambda_infty.png"),
       plot = myplt, width = 12, height = 4, dpi = 600)
saveRDS(list(seed = 1L, Ns = Ns, noise = noise, circle = results1,
             interval = results2, sphere = results3),
        file.path(figure_dir, "dist_lambda_infty_results.rds"))
