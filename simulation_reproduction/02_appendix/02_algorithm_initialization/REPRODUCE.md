# Reproduce the initialization figures

- `iso_1d_2d.png`: run `code/isomap_initialization.R`. The left panel shows the 1D ISOMAP initialization in dark green; the right panel shows the 2D ISOMAP initialization in dashed blue. Both panels overlay the estimate on the lower-opacity dark-gray true manifold and intentionally omit panel labels so these can be supplied in the figure caption.
- `isomap_unwind.png` and `isomap_unwind_first_vs_converged.png`: run `code/isomap_unwinding_source.R`. This standalone script is extracted from the historical unwinding section and uses the same fixed seed and cashew-shaped sample for both figures. In the comparison figure, the left panel is the initial fit and the right panel is the converged fit; solid red denotes Naive initialization and dashed blue denotes ISOMAP initialization.
- `example_star_convex.png`: this is an Ipe drawing rather than a simulation output; retain the supplied manuscript asset. It has no R generator.

Install CompactPME, `vegan`, and `scico` before running the ISOMAP simulation scripts.
