library(pracma)
library(dplyr)
library(ggplot2)
library(TDA)
library(plot3D)
library(tidyr)
library(patchwork)
library(png)
library(grid)
library(gridExtra)

# 2d TDA
set.seed(1)
fibonaccisphere = function(n=1000, r=1, out.xyz=TRUE, out.sph=FALSE) {
  
  if (n<1 | round(n)!=n) stop('n must be a positive integer')
  if (!out.xyz & !out.sph) stop('either out.xyz and/or out.sph must be TRUE')
  
  goldenratio = (1+sqrt(5))/2
  i = seq(n)-0.5
  z = 1-2*i/n # z-coordinate for unit sphere
  theta = acos(pmax(-1,pmin(1,z))) # polar angle
  phi = (2*pi*i/goldenratio)%%(2*pi)
  
  if (out.xyz) {
    x = r*sin(theta)*cos(phi)
    y = r*sin(theta)*sin(phi)
    z = r*z
    if (out.sph) {
      out = cbind(x=x,y=y,z=z,theta=theta,phi=phi)
    } else {
      out = cbind(x=x,y=y,z=z)
    }
  } else {
    out = cbind(theta=theta,phi=phi)
  }
  
  return(out)
  
}

runTDA_2d3D <- function(data, alpha=0.01, boot_size=30)
{
  Xlim <- c(min(data[,1]), max(data[,1]))
  Ylim <- c(min(data[,2]), max(data[,2]))
  Zlim <- c(min(data[,3]), max(data[,3]))
  lim <- cbind(Xlim, Ylim,Zlim)
  by <- 0.05
  
  Diag <- gridDiag(data, distFct, lim = lim, by = by)
  boot <- bootstrapDiagram(X=data, FUN=distFct, lim = lim,by=by, parallel=T,alpha=alpha,B=boot_size)
  
  list(Diag=Diag, boot=boot, alpha =alpha)
}

tda_2_plot2d <- function(tda_res, label="", hide_legend=F, title=""){
  Diag <- tda_res$Diag
  boot <- tda_res$boot
  alpha <- tda_res$alpha
  
  conf_label <- paste0(round((1 - alpha) * 100), "% confidence band")
  
  pd_df <- Diag$diagram %>% unclass() %>% as.data.frame()
  
  grid_points <- seq(min(pd_df$Birth)-1, max(pd_df$Birth)+5, length.out = 100)
  
  #ribbon_df <- data.frame(
  #  x = c(grid_points - boot, rev(grid_points + boot)),
  #  y = c(grid_points + boot, rev(grid_points - boot)),
  #  Band=factor(conf_label))
  
  ribbon_df <- data.frame(
    x= c(grid_points, rev(grid_points)),
    y=c(grid_points + 2*boot,rev(grid_points)),
    Band = conf_label
  )
  
  
  xmax <- max(pd_df$Birth)*1.05
  ymax <- max(pd_df$Death)*1.05
  
  legendpos <- ifelse(hide_legend, "none", "bottom")
  
  pd_df %>%  
    ggplot() + 
    geom_point(aes(x=Birth,y=Death, color = factor(dimension), shape=factor(dimension)),size=2) +
    scale_color_manual(name=paste0(label, "Dimension"),values=c("0"="black", "1"="red", "2"= "blue")) +
    scale_shape_manual(name=paste0(label, "Dimension"),values = c("0" = 16, "1" = 2, "2"= 5)) +
    geom_abline(slope=1,intercept=0)  + 
    geom_polygon(data=ribbon_df, aes(x=x,y=y, fill = Band),alpha=0.2) + 
    theme_bw() + 
    theme(panel.grid=element_blank(), legend.position=legendpos,legend.background=element_rect(fill="white",color="black"))+ 
    scale_fill_manual(values = setNames("red", levels(ribbon_df$Band)), name = paste0(label, "Band"))+theme(aspect.ratio=1)+
    #xlim(c(0,xmax))+ylim(c(0,ymax))
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) + ggtitle(title)
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
# Load / setup data
tda_sphere <- readRDS("results_sph/2d_sphere_tda_half.rds")
(p_tda_sph <- tda_2_plot2d(tda_sphere,title="(iv)",hide_legend=T))
#set.seed(1)
ggsave("panel_figs/panel_A_4.png", width=3200, height=3200, units = "px",dpi=600)







# 2d Viz plots

genSurface <- function(f){
  M <- 50
  thetas <- seq(from = 0, to = 2*pi, length.out=M)
  phis <- seq(from =-pi/2, to = pi/2, length.out=M)
  rr <- 1
  grid <- expand.grid(thetas,phis,1)
  cartgrid <- t(apply(grid, 1, sph2cart))
  
  # Evaluate manifold
  fvalues <- f(cartgrid)
  
  Xmat <- matrix(fvalues[,1], nrow = M, ncol = M, byrow = TRUE)
  Ymat <- matrix(fvalues[,2], nrow = M, ncol = M, byrow = TRUE)
  Zmat <- matrix(fvalues[,3], nrow = M, ncol = M, byrow = TRUE)
  
  list(X=Xmat,Y=Ymat,Z=Zmat)
}

plot_datafit_2d3D_static <- function(X, fit, theta=40,phi=20, bty="n", title="Data and Surface",
                                     point_size = 0.5){
  surf <- genSurface(fit)
  
  
  
  # Get the lim given by min over all dimensions and max over all dimensions
  #lims <- c(min(extract_lims[1,]), max(extract_lims[2,]))
  lims <- c(min(X), max(X))
  
  # Start with data points
  scatter3D(
    X[,1], X[,2], X[,3],
    pch = 16, cex = point_size, col = "black",
    theta = theta, phi = phi,
    main = title,
    bty = bty,
    ticktype = "detailed", xlim=lims, ylim=lims,zlim=lims
    
    
  )
  
  # Overlay surface
  surf3D(
    surf$X, surf$Y, surf$Z,
    add = TRUE,
    colvar = surf$Z, # overlay on existing plot
    # col = "Viridis",
    alpha = 0.6,             # semi-transparent
    border = "grey",
    #xlim=c(-5,5),
    bty = bty,
    colkey=F,
    
    xlim=lims, ylim=lims,zlim=lims
    
  )
  
  
}

pme_sphere <- readRDS("results_sph/2d_sphere_pme_half.rds")

# choose lambda



# get plots of three choices of lambda
get_closest_lambda_ind <- function(result, lambda_target) {
  
  idx <- sapply(lambda_target, function(l) {
    which.min(abs(result$lambda_list - l))
  })
  
  unique(idx)
}




# Combine TDA and Viz






## CV Plots

#stackoverflow
# Source - https://stackoverflow.com/a/34205938
# Posted by stas g
# Retrieved 2026-02-25, License - CC BY-SA 3.0



# Plot all kernels or just pelletier?
cv_plots2d <- function(result, f_true, cutoff= 5e-4){
  
  
  reshapedf <- function(df){
    as.data.frame(df)  %>% pivot_longer(cols=-lambda,names_to="kernel_method",values_to="value")
  }

    # Results for pelletier only
    
    meanplt <- result$mean %>% reshapedf() %>% filter(kernel_method=="pelletier") %>% ggplot(aes(x=lambda, y=value)) + geom_point()+ geom_line()+ scale_x_log10() + theme_bw()+
      theme(panel.grid=element_blank()) + ylab("Mean Squared Residual") + xlab("λ")+ ggtitle("(ii)") + geom_vline(xintercept=cutoff, linetype="dashed")
    
    
    sdplt <- result$sd %>% reshapedf() %>% filter(kernel_method=="pelletier") %>% ggplot(aes(x=lambda, y=value)) + geom_point()+ geom_line()+ scale_x_log10() + theme_bw()+
      theme(panel.grid=element_blank()) + ylab("Standard Deviation") + xlab("λ")+ ggtitle("(iii)") 
    
    # Min_df
    mindf <- result$cv %>% reshapedf() %>% filter(kernel_method=="pelletier") %>% group_by(kernel_method) %>% slice_min(value, with_ties=F)
    # thresholded min
    mindf2 <- result$cv %>% reshapedf() %>% filter(lambda <= cutoff, kernel_method=="pelletier") %>% group_by(kernel_method) %>% slice_min(value, with_ties=F)
    
    mindfs <- rbind(cbind(mindf, optimal="naive"), cbind(mindf2, optimal="thresholded"))
    
    
    cvplt <-result$cv %>% reshapedf() %>%  filter(kernel_method=="pelletier", lambda <= cutoff) %>% 
      ggplot(aes(x=lambda, y=value)) + geom_point()+ geom_line()+ scale_x_log10() + 
      theme_bw()+ geom_vline(xintercept = mindf2$lambda, col="red",linetype="dashed")+
      theme(panel.grid=element_blank(), legend.position="none")+
      ylab("Coefficient of Variation") + xlab("λ") + ggtitle("(i)")
    
    
  
  #(cvplt + meanplt)/ (sdplt+ plot_spacer()) 
  list(cv = cvplt, mean = meanplt, sd = sdplt)
}

ftrue <- flower2d3D_unif(1000)
cv_plots2d(pme_sphere, ftrue$X)
#ggsave("panel_figs/cv_res.png", width=4000, height=3200, units = "px",dpi=300)
cvplots <- cv_plots2d(pme_sphere, ftrue$X)

cvplots$cv
ggsave("panel_figs/panel_A_1.png", width=3200, height=3200, units = "px",dpi=600)
cvplots$mean
ggsave("panel_figs/panel_A_2.png", width=3200, height=3200, units = "px",dpi=600)
cvplots$sd
ggsave("panel_figs/panel_A_3.png", width=3200, height=3200, units = "px",dpi=600)


# Plots of surfaces with optimal lambdas ======
get_closest_lambda_ind <- function(result, lambda_target) {
  
  idx <- sapply(lambda_target, function(l) {
    which.min(abs(result$lambda_list - l))
  })
  
  idx
}

reshapedf <- function(df){
  as.data.frame(df)  %>% pivot_longer(cols=-lambda,names_to="kernel_method",values_to="value")
}

# Min_df

cutoff <- 5e-4
# Min_df
mindf <- pme_sphere$cv %>% reshapedf() %>% filter(kernel_method=="pelletier") %>% group_by(kernel_method) %>% slice_min(value, with_ties=F)
# thresholded min
mindf2 <- pme_sphere$cv %>% reshapedf() %>% filter(lambda <= cutoff, kernel_method=="pelletier") %>% group_by(kernel_method) %>% slice_min(value, with_ties=F)

mindfs <- rbind(cbind(mindf, type="naive"), cbind(mindf2, type="cutoff"))

lam_min <- mindfs$lambda
kern_min <- mindfs$type

lam_min_ind <- get_closest_lambda_ind(pme_sphere, lam_min)

ilabs <- c("(ii)", "(iii)")
for(i in 1:length(lam_min)){
  
  png(paste0("panel_figs/panel_B_",i,".png"), width = 3200, height=3200, res=600)
  par(mar=c(0.02,0.1,3,0.1))
  
  lam_i <- lam_min_ind[i]
  mytitle <- paste0(ilabs[i]," ", kern_min[i], ", λ = ",formatC(pme_sphere$lambda_list[lam_i],
                                                                format = "e", digits=2))
  plot_datafit_2d3D_static(pme_sphere$dataset, pme_sphere$spline_list[[lam_i]], phi=40, point_size=0.4,
                           title=mytitle)
  
  dev.off()
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






lams <- c(10^-12, mindf2$lambda ,10^-2)
lams_ind <- get_closest_lambda_ind(pme_sphere, lams)

plotz <- vector("list", length=length(lams))

prelab <- c("(vi)", "(vii)", "(viii)")

for(i in 1:length(lams)){
  
  png(paste0("panel_figs/panel_B_",i+1,".png"), width = 3200, height=3200, res=600)
  par(mar=c(0.02,0.1,3,0.1))
  lam_i <- lams_ind[i]
  plot_datafit_2d3D_static(pme_sphere$dataset, pme_sphere$spline_list[[lam_i]], phi=40, point_size=0.4,
                           title=paste0(prelab[i], " λ = ",formatC(pme_sphere$lambda_list[lam_i],
                                                            format = "e", digits=2)))
  
  dev.off()
}



png(paste0("panel_figs/panel_B_1.png"), width = 3200, height=3200, res=600)
ftrue <- function(t) flower2d3D_func(t, petals=5)
par(mar=c(0.02,0.1,3,0.1))
plot_datafit_2d3D_static(pme_sphere$dataset, ftrue, phi=40, point_size=0.4, title="(v) Data and True Surface")
dev.off()
#plot_datafit_2d3D(pme_sphere$dataset, ftrue)



filesA <- sprintf("panel_figs/panel_A_%d.png", 1:4)
filesB <- sprintf("panel_figs/panel_B_%d.png", 1:4)
files <- c(filesA, filesB)
# Read all images
imgs <- lapply(files, function(f) {
  rasterGrob(readPNG(f), interpolate = TRUE)
})
# Arrange in a 2x2 grid
png("panel_figs/panel8_fig.png", width = 12800, height = 6400,res=600)
grid.arrange(grobs = imgs, ncol = 4, nrow = 2)
dev.off()








