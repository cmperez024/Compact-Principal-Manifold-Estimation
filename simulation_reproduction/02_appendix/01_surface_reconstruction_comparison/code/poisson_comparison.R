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
# Compute this study and retain its original manuscript surface rendering.
options(rgl.useNULL = TRUE)
required <- c("SurfaceReconstruction", "TDA", "FNN", "rgl", "png")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " first.")
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
  k <- 10
  d <- FNN::get.knn(pts, k)$nn.dist[,k]       # k-th NN distance
  w <- d^3                               # local volume ~ r^3
  prob <- w / max(w)
  
  keep <- sample(seq_len(nrow(pts)), N, prob=prob)
  
  list(X = pts[keep,] + rnorm(3*N, sd=noise),
       t = S[keep,])
}


set.seed(1)
n <- if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) 100L else 1000L
X <- unifmoon2d3D(n, bend = 1.5*pi, bend_radius = 2, noise = 0.05)$X
limits <- apply(X, 2L, range)
by <- if (n == 100L) 0.3 else 0.1
tda <- TDA::gridDiag(X = X, FUN = TDA::distFct, lim = limits, by = by)
png(file.path(figure_dir, "mytda.png"), width = 1200, height = 1200, res = 200)
plot(tda$diagram, asp = 1)
legend(1.5, 0.75, c("0", "1", "2"), col = c(1, 2, 4),
       pch = c(20, 2, 5), cex = 0.9, pt.lwd = 2)
dev.off()

mesh <- SurfaceReconstruction::PoissonReconstruction(X)
# The final manual camera rotation was not recorded in the notebook. Keep
# this workflow's original shaded snapshot for the fixed manuscript study.
# The freshly computed mesh is saved separately rather than pretending that
# an approximate camera or a changed CGAL triangulation matches the original.
render_mode <- Sys.getenv("PME_POISSON_RENDER", unset =
  if (identical(Sys.getenv("PME_SMOKE_TEST"), "1")) "rgl" else "manuscript")
if (!render_mode %in% c("manuscript", "rgl")) {
  stop("PME_POISSON_RENDER must be 'manuscript' or 'rgl'.")
}
surface_image <- file.path(figure_dir, "myhole.png")
if (render_mode == "manuscript") {
  reference_image <- file.path(script_dir, "..", "assets",
                               "poisson_surface_manuscript.png")
  if (!file.exists(reference_image)) {
    stop("Missing manuscript surface asset: ", reference_image,
         ". Set PME_POISSON_RENDER=rgl for a fresh shaded reconstruction.")
  }
  if (!file.copy(reference_image, surface_image, overwrite = TRUE)) {
    stop("Could not copy the manuscript surface snapshot.")
  }
} else {
  if (!requireNamespace("webshot2", quietly = TRUE)) {
    stop("Install webshot2 for automated rgl rendering.")
  }
  device <- rgl::open3d(windowRect = c(0, 0, 1000, 1000))
  rgl::bg3d(color = "white")
  rgl::shade3d(mesh, color = "pink", lit = TRUE, smooth = TRUE,
               specular = "white", shininess = 50)
  rgl::wire3d(mesh, color = "black")
  # Upright U-shaped orientation. This is explicitly a fresh rendering;
  # the exact original camera is preserved only by the manuscript snapshot.
  rgl::par3d(userMatrix = rgl::rotationMatrix(-pi/2, 0, 0, 1),
             FOV = 30, zoom = 0.9)
  rgl::snapshot3d(surface_image, width = 1000, height = 1000, webshot = TRUE)
  rgl::close3d(device)
  if (!file.exists(surface_image) || file.info(surface_image)$size == 0) {
    stop("The shaded rgl snapshot failed.")
  }
}

img1 <- png::readPNG(file.path(figure_dir, "mytda.png"))
img2 <- png::readPNG(file.path(figure_dir, "myhole.png"))
png(file.path(figure_dir, "combined.png"), width = 2200, height = 1200)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
grid::grid.raster(img1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
grid::grid.raster(img2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
saveRDS(list(seed = 1L, data = X, diagram = tda$diagram, mesh = mesh,
             surface_image_mode = render_mode),
        file.path(figure_dir, "poisson_comparison_results.rds"))
