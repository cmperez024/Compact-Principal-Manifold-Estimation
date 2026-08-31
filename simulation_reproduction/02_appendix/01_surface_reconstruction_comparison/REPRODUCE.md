# Reproduce `combined.png`

Run `code/poisson_comparison.Rmd` in order. It generates the point cloud, computes its persistence diagram, reconstructs the Poisson surface, and combines the two panels.

The notebook expects `rgl`, `plotly`, `FNN`, `TDA`, `png`, `grid`, and the
`SurfaceReconstruction` package. Install the latter from its upstream GitHub
repository before running the notebook. The historical `rgl.snapshot` step was
interactive; execute it after positioning the mesh view to create `myhole.png`
before the combination chunk.
