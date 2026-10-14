# Reproduce `combined.png`

Run `code/poisson_comparison.R`. It generates the point cloud, computes its
persistence diagram, reconstructs the Poisson surface, and combines the panels.

The script requires `FNN`, `TDA`, `png`, `rgl`, and `SurfaceReconstruction`
(from `stla/SurfaceReconstruction`). The default publication run retains the
original pink, shaded surface with black mesh edges using this workflow's own
`assets/poisson_surface_manuscript.png`. The final manual camera rotation was
not recorded in the historical notebook; the original snapshot preserves it
exactly. The newly computed mesh is saved in `poisson_comparison_results.rds`
and is not misrepresented as the mesh underlying the historical snapshot.

The script writes `mytda.png`, `myhole.png`, and `combined.png`. Its default
combined output has been checked pixel-for-pixel against the original
manuscript figure. The persistence device is closed before reading the image.
No other simulation script or another simulation's results are used.

For a fresh shaded rgl rendering, set `PME_POISSON_RENDER=rgl`. This requires
`webshot2` and Chrome/Chromium. It uses pink smooth shading, white specular
highlights, black mesh edges, and an upright U-shaped view. Its camera and
CGAL triangulation may differ from the original snapshot. Reduced tests also
use this mode so a fresh mesh is shown for the reduced sample.

Paths resolve from the script. `PME_OUTPUT_DIR` selects an alternate output
directory and `PME_SMOKE_TEST=1` selects a reduced validation sample.

With current RcppCGAL headers, building SurfaceReconstruction requires C++17.
For Windows source installation, set `CXX14STD = -std=gnu++17` in a temporary
Makevars file and point `R_MAKEVARS_USER` to that file during installation.
