# Reproduce the initialization figures

- `iso_1d_2d.png`: run the **Circle** chunk in `code/isomap_initialization.Rmd`.
- `isomap_unwind_first_vs_converged.png`: run **Isomap unwound example** through its `ggsave` call in `code/isomap_unwinding_source.Rmd`.
- `isomap_unwind.png` and `example_star_convex.png`: no standalone generator with those filenames exists in the workspace. They are manuscript-source illustrations; the former corresponds to the original-data/ISOMAP diagnostic in the unwinding section.

Load CompactPME plus `vegan` and `ggplot2` before running these sections.
