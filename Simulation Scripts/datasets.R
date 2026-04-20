# Fibonacci, rotation2, and rotation 3 from cooltools package
fibonaccisphere = function(n=1000, r=1) {
  
  if (n<1 | round(n)!=n) stop('n must be a positive integer')
  
  goldenratio = (1+sqrt(5))/2
  i = seq(n)-0.5
  z = 1-2*i/n # z-coordinate for unit sphere
  phi = acos(pmax(-1,pmin(1,z))) # polar angle
  theta = (2*pi*i/goldenratio)%%(2*pi)
  
  x = r*sin(phi)*cos(theta)
  y = r*sin(phi)*sin(theta)
  z = r*z
  out = cbind(x=x,y=y,z=z)
  
  return(out)
  
}


rotation2 = function(angle) {
  c = cos(angle)
  s = sin(angle)
  return(rbind(c(c,-s),c(s,c)))
  
}

rotation3 = function(u, angle=NULL) {
  n = length(u)
  if (n==3) {
    unorm = sqrt(sum(u^2))
    if (unorm==0) {
      R = diag(3)
    } else {
      if (is.null(angle)) angle = unorm
      u = u/unorm
      c = cos(angle)
      s = sin(angle)
      R = rbind(c(c+u[1]^2*(1-c), u[1]*u[2]*(1-c)-u[3]*s, u[1]*u[3]*(1-c)+u[2]*s),
                c(u[2]*u[1]*(1-c)+u[3]*s, c+u[2]^2*(1-c), u[2]*u[3]*(1-c)-u[1]*s),
                c(u[3]*u[1]*(1-c)-u[2]*s, u[3]*u[2]*(1-c)+u[1]*s, c+u[3]^2*(1-c)))
    }
    
    return(R)
  } else {
    stop('u must be a 3-element vector')
  }
}


moon1d2D <- function(N, inner_radius=1, outer_radius=2, noise=0)
{
  # Parameters (example)
  r1 <- inner_radius    # inner radius
  r2 <- outer_radius      # outer radius (must be > r1)
  stopifnot(r2 > r1)
  #n <- N/4 # points per segment (increase for smoother curve)
  nmain <- round(4*N/5)
  ncap <- round(N/5)
  
  
  r_cap <- (r2 - r1) / 2
  k <- (r1 + r2) / 2    # midpoint y for cap centers
  
  # 1) Outer left semicircle: top -> bottom (theta from pi/2 to 3pi/2)
  theta_outer <- seq(pi/2, 3*pi/2, length.out = nmain)
  x_outer <- r2 * cos(theta_outer)
  y_outer <- r2 * sin(theta_outer)
  
  # 2) Bottom right-half cap: connects (0, -r2) -> (0, -r1)
  phi_bot <- seq(-pi/2, pi/2, length.out = ncap)           # right half (x >= 0)
  x_bot_cap <- r_cap * cos(phi_bot)
  y_bot_cap <- -k + r_cap * sin(phi_bot)                # center at (0, -k)
  
  # 3) Inner left semicircle: bottom -> top (theta from 3pi/2 -> pi/2)
  theta_inner <- seq(3*pi/2, pi/2, length.out = nmain)
  x_inner <- r1 * cos(theta_inner)
  y_inner <- r1 * sin(theta_inner)
  
  # 4) Top right-half cap: connects (0, r1) -> (0, r2)
  phi_top <- seq(-pi/2, pi/2, length.out = ncap)           # right half (x >= 0)
  x_top_cap <- r_cap * cos(phi_top)
  y_top_cap <-  k + r_cap * sin(phi_top)                # center at (0, k)
  
  # Combine in order to form the closed loop
  x <- c(x_outer, x_bot_cap, x_inner, x_top_cap)
  y <- c(y_outer, y_bot_cap, y_inner, y_top_cap)
  
  cbind(x,y) + rnorm(2*length(x), sd=noise)
}

moon1d3D <- function(N, inner_radius=1, outer_radius=2, noise=0)
{
  # Parameters (example)
  r1 <- inner_radius    # inner radius
  r2 <- outer_radius      # outer radius (must be > r1)
  stopifnot(r2 > r1)
  #n <- N/4 # points per segment (increase for smoother curve)
  nmain <- round(4*N/5)
  ncap <- round(N/5)
  
  
  r_cap <- (r2 - r1) / 2
  k <- (r1 + r2) / 2    # midpoint y for cap centers
  
  # 1) Outer left semicircle: top -> bottom (theta from pi/2 to 3pi/2)
  theta_outer <- seq(pi/2, 3*pi/2, length.out = nmain)
  x_outer <- r2 * cos(theta_outer)
  y_outer <- r2 * sin(theta_outer)
  
  # 2) Bottom right-half cap: connects (0, -r2) -> (0, -r1)
  phi_bot <- seq(-pi/2, pi/2, length.out = ncap)           # right half (x >= 0)
  x_bot_cap <- r_cap * cos(phi_bot)
  y_bot_cap <- -k + r_cap * sin(phi_bot)                # center at (0, -k)
  
  # 3) Inner left semicircle: bottom -> top (theta from 3pi/2 -> pi/2)
  theta_inner <- seq(3*pi/2, pi/2, length.out = nmain)
  x_inner <- r1 * cos(theta_inner)
  y_inner <- r1 * sin(theta_inner)
  
  # 4) Top right-half cap: connects (0, r1) -> (0, r2)
  phi_top <- seq(-pi/2, pi/2, length.out = ncap)           # right half (x >= 0)
  x_top_cap <- r_cap * cos(phi_top)
  y_top_cap <-  k + r_cap * sin(phi_top)                # center at (0, k)
  
  # Combine in order to form the closed loop
  x <- c(x_outer, x_bot_cap, x_inner, x_top_cap)
  y <- c(y_outer, y_bot_cap, y_inner, y_top_cap)
  
  t0 <- seq(0,1,length.out=length(x))
  cbind(x,y, cos(4*pi*t0)) + rnorm(3*length(x), sd=noise)
}

# 1d curve in 2D
flower1d2D <- function(n, petals, noise, t){
  t0 <- if(missing(t)) sort(runif(n)) else t
  
  r <- (1 + 0.3 * sin(petals *2*pi* t0))
  return(list(t=t0, X = cbind(r * cos(2*pi*t0),  r * sin(2*pi*t0))+rnorm(2*n, sd=noise)))
}

halfflower1d2D <- function(n, petals, noise, t, rotation=0){
  t0 <- if(missing(t)) sort(runif(n)) else t
  
  r <- (1 + 0.3 * sin(petals *pi* t0))
  return(list(t=t0, X = cbind(r * cos(pi*t0),  r * sin(pi*t0)) %*% rotation2(rotation) +rnorm(2*n, sd=noise)))
} 


flower1d3D <- function(n, petals, noise, t){
  t0 <- if(missing(t)) sort(runif(n)) else t
  
  r <- (1 + 0.3 * sin(petals *2*pi* t0))
  return(list(t=t0, X = cbind(r * cos(2*pi*t0),  r * sin(2*pi*t0), sin(3*pi*t0))+rnorm(3*n, sd=noise)))
}


circle1d2D <- function(n, xscale=1, yscale=1, center=c(0,0), noise=0, t){
  t0 <- if(missing(t)) sort(runif(n)) else t
  #print(str(t0))
  
  X <- cbind(center[1] + xscale*cos(2*pi*t0), center[2]+yscale*sin(2*pi*t0)) + rnorm(2*n, sd=noise)
  return(list(t=t0,X=X))
  
}

subcircle1d2D <- function(n, xscale=1, yscale=1, center=c(0,0), noise=0, t){
  t0 <- if(missing(t)) sort(runif(n)) else t
  #print(str(t0))
  
  X <- cbind(center[1] + xscale*cos(1.5*pi*t0), center[2]+yscale*sin(1.5*pi*t0)) + rnorm(2*n, sd=noise)
  return(list(t=t0,X=X))
  
}


sine1d2D <- function(n, rotation=0, period=2*pi, noise=0, tval){
  t0 <- if(missing(tval)) sort(runif(n)) else tval
  
  list(X=cbind(t0, sin(period*t0)) %*% t(rotation2(rotation)) + rnorm(2*n, sd=noise),
       t = t0)
}

random_sphere <- function(n) {
  u <- runif(n, -1, 1)
  phi <- runif(n, 0, 2*pi)
  theta <- acos(u)     # correct relation
  x <- sin(theta)*cos(phi)
  y <- sin(theta)*sin(phi)
  z <- cos(theta)
  cbind(x,y,z)
}

spiral1d2D <- function(n, period = 4*pi, scale = 1, noise=0){
  
  t0 <- sort(runif(n))
  
  list(X=cbind(scale*t0*cos(period*t0), scale*t0*sin(period*t0))+ rnorm(2*n, sd=noise),
       t = t0)
}


#chatgpt
spiral1d2D_arclength <- function(n, period = 4*pi, scale = 1, noise = 0){
  
  omega <- period
  
  # arc length function
  s_fun <- function(t){
    0.5 * t * sqrt(1 + omega^2 * t^2) +
      asinh(omega * t) / (2 * omega)
  }
  
  # total length
  L <- s_fun(1)
  
  # sample uniformly in arc length
  s_star <- sort(runif(n, 0, L))
  
  # invert numerically
  t_star <- sapply(s_star, function(s_target){
    uniroot(function(t) s_fun(t) - s_target,
            interval = c(0,1))$root
  })
  
  # evaluate curve
  X <- cbind(scale * t_star * cos(omega * t_star),
             scale * t_star * sin(omega * t_star))
  
  # add noise
  if(noise > 0){
    X <- X + matrix(rnorm(2*n, sd=noise), ncol=2)
  }
  
  list(X = X, t = t_star)
}


sine1d3D <- function(n, rotation=0, period=2*pi, noise=0, tval){
  t0 <- if(missing(tval)) sort(runif(n)) else tval
  X <- t(cbind(t0, sin(period*t0),cos(period*t0)))
  X <- rotation3(c(1,1,0.2),rotation)%*% X
  X <- t(X) + rnorm(3*n, sd=noise)
  list(X=X,t=t0)
}

spiral1d2D <- function(n, period = 4*pi, scale = 1, noise=0){
  
  t0 <- sort(runif(n))
  
  list(X=cbind(scale*t0*cos(period*t0), scale*t0*sin(period*t0))+ rnorm(2*n, sd=noise),
       t = t0)
}

spiral1d3D <- function(n, period = 4*pi, scale = 1, noise=0, rotation=0){
  
  t0 <- sort(runif(n))
  
  list(X=cbind(scale*t0*cos(period*t0), scale*t0*sin(period*t0), scale*t0)+ rnorm(3*n, sd=noise),
       t = t0)
}

helix1d3D <- function(n, period = 4*pi, xyscale = 1, zscale=1, noise=0, rotation=0){
  
  t0 <- sort(runif(n))
  
  list(X=cbind(xyscale*cos(period*t0), xyscale*sin(period*t0), zscale*t0) %*% t( rotation3(c(1,1,0.2),rotation))+ rnorm(3*n, sd=noise),
       t = t0)
}



# axisscale is c(xr,yr,zr) where xr multiplies x axis, etc
ellipse2d3D <- function(n, axisscale, noise){
  # scale x by xr, y by yr, z by zr 
  X <- sweep(fibonaccisphere(n), 2, axisscale, "*")
  X + rnorm(3*n, sd=noise)
}

unifmoon2d3D <- function(N, oversample = 5, bend = pi,
                         radius = 1, bend_radius = 1.2, noise=0) {
  
  # 1. oversample sphere uniformly
  S <- random_sphere(N * oversample)
  
  x <- radius * S[,1]
  y <- radius * S[,2]
  z <- radius * S[,3]
  
  # 2. bend map
  new_x <- cos(bend * x / (2*radius)) * (bend_radius * radius + y)
  new_y <- sin(bend * x / (2*radius)) * (bend_radius * radius + y)
  new_z <- z
  
  pts <- cbind(new_x, new_y, new_z)
  
  # 3. Rejection thinning for uniform coverage
  # compute kernel density to approximate sampling density
  library(FNN)
  k <- 10
  d <- get.knn(pts, k)$nn.dist[,k]       # k-th NN distance
  w <- d^3                               # local volume ~ r^3
  prob <- w / max(w)
  
  keep <- sample(seq_len(nrow(pts)), N, prob=prob)
  
  list(X = pts[keep,] + rnorm(3*N, sd=noise),
       t = S[keep,])
}


# 
# macaraoni / moon data =======
moon2d3D <- function(n_theta = 80, n_phi = 80, bend = 2*pi, radius = 1, bend_radius = 1.2, noise=0) {
  
  # equal-area theta sampling (reduces pole clustering)
  u <- runif(n_theta, -1, 1)
  theta <- asin(u)                # maps -1..1 -> -pi/2..pi/2
  phi <- seq(0, 2*pi, length.out = n_phi)
  
  #theta <- seq(-pi/2, pi/2, length.out = n_theta)
  # phi <- seq(0, 2*pi, length.out = n_phi)
  g <- expand.grid(theta = theta, phi = phi)
  
  x <- radius * cos(g$theta) * cos(g$phi)
  y <- radius * cos(g$theta) * sin(g$phi)
  z <- radius * sin(g$theta)
  
  new_x <- cos(bend * x / (2*radius)) * (bend_radius * radius + y)
  new_y <- sin(bend * x / (2*radius)) * (bend_radius * radius + y)
  new_z <- z
  
  
  
  pts <- cbind(x= new_x, y=new_y, z=new_z)+ rnorm(3*length(new_x), sd=noise)
  
  #pts
  g <- as.matrix(g)
  return(list(t = sph2cart(cbind(g,1)), X=pts))
}



# chatgpt. less sampling in poles
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

flower2d3D_unif <- function(N = 300,
                            r0 = 1, a = 0.3, petals = 6,
                            b = 0.5, noise = 0.01,
                            t = NULL) {
  
  # If no sphere provided, generate uniform Fibonacci sphere
  if (is.null(t)) {
    t <- fibonaccisphere(N)  # N x 3 matrix, uniform on S^2
  }
  
  # Ensure unit sphere
  t <- t / sqrt(rowSums(t^2))
  
  # Convert to spherical coordinates
  coords <- cart2sph(t)
  phi   <- coords[, "theta"]     # longitude
  elev  <- coords[, "phi"]       # elevation
  z     <- t[, 3]                # vertical coordinate
  
  # Radial modulation in xy-plane
  r_xy <- r0 * (1 + a * cos(petals * phi))
  
  # Flower deformation
  x <- r_xy * sqrt(1 - z^2) * cos(phi)
  y <- r_xy * sqrt(1 - z^2) * sin(phi)
  z <- b * z
  
  pts <- cbind(x, y, z) + matrix(rnorm(3*N, sd = noise), N, 3)
  
  return(list(t = t, X = pts))
}

flower2d3D_func <- function(t,
                            r0 = 1, a = 0.3, petals = 6,
                            b = 0.5
) {
  
  # Ensure unit sphere
  t <- t / sqrt(rowSums(t^2))
  
  # Convert to spherical coordinates
  coords <- cart2sph(t)
  phi   <- coords[, "theta"]     # longitude
  elev  <- coords[, "phi"]       # elevation
  z     <- t[, 3]                # vertical coordinate
  
  # Radial modulation in xy-plane
  r_xy <- r0 * (1 + a * cos(petals * phi))
  
  # Flower deformation
  x <- r_xy * sqrt(1 - z^2) * cos(phi)
  y <- r_xy * sqrt(1 - z^2) * sin(phi)
  z <- b * z
  
  cbind(x, y, z)
}