# Reproduce the four-by-five density-and-curve figure for the
# quantile-calibrated numerical-nonidentifiability example in Section 5.1.
# The implementation uses only base R.

script_directory <- function() {
  file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_argument)) {
    return(dirname(normalizePath(sub("^--file=", "", file_argument[[1L]]))))
  }
  normalizePath(getwd())
}

output_file <- Sys.getenv(
  "PME_FIGURE_OUTPUT",
  unset = file.path(
    dirname(script_directory()),
    "figure_5_1_numerical_nonidentifiability.png"
  )
)

sigma_circle <- 0.1
plot_limit <- 1.5
grid_size <- 501L
radial_quadrature_size <- 4096L
density_ceiling <- 0.65

rice_quantile <- function(probability) {
  sigma_circle * sqrt(
    qchisq(
      probability,
      df = 2,
      ncp = 1 / sigma_circle^2
    )
  )
}

curve_parameters <- function(q) {
  eta <- 1 / (10 * (q + 1))
  lower_radius <- rice_quantile(eta)
  upper_radius <- rice_quantile(1 - eta)

  list(
    q = q,
    eta = eta,
    lower_radius = lower_radius,
    upper_radius = upper_radius,
    radial_width = upper_radius - lower_radius,
    kappa = 10L * q,
    sigma_q = sigma_circle / sqrt(q)
  )
}

# Add possibly repeated weighted indices to a dense vector. rowsum performs
# the aggregation in compiled code and is substantially faster than an R loop.
add_mass <- function(target, index, weight) {
  aggregated <- rowsum(weight, index, reorder = FALSE)
  occupied <- as.integer(rownames(aggregated))
  target[occupied] <- target[occupied] + aggregated[, 1L]
  target
}

# Deposit a deterministic approximation to the curve-supported latent
# distribution on
# the plotting grid. The radial midpoint rule is applied in Rice-probability
# coordinates, so every generated source point has the same probability mass.
deposit_curve_distribution <- function(parameters, grid, quadrature_size) {
  n_grid <- length(grid)
  spacing <- grid[2L] - grid[1L]
  lower_grid_limit <- grid[1L]
  kappa <- parameters$kappa

  probability <- parameters$eta +
    (1 - 2 * parameters$eta) *
    ((seq_len(quadrature_size) - 0.5) / quadrature_size)
  radius <- rice_quantile(probability)
  phase_argument <- 1 -
    2 * (radius - parameters$lower_radius) / parameters$radial_width
  phase <- acos(pmax(-1, pmin(1, phase_argument)))

  cell <- 0:(kappa - 1L)
  theta_plus <- outer(
    phase / kappa,
    2 * pi * cell / kappa,
    "+"
  )
  theta_minus <- outer(
    -phase / kappa,
    2 * pi * (cell + 1L) / kappa,
    "+"
  )
  theta <- c(theta_plus, theta_minus)
  radius_repeated <- rep(radius, times = 2L * kappa)

  source_x <- radius_repeated * cos(theta)
  source_y <- radius_repeated * sin(theta)
  grid_x <- (source_x - lower_grid_limit) / spacing + 1
  grid_y <- (source_y - lower_grid_limit) / spacing + 1
  ix <- floor(grid_x)
  iy <- floor(grid_y)
  fraction_x <- grid_x - ix
  fraction_y <- grid_y - iy

  stopifnot(
    all(ix >= 1L), all(ix < n_grid),
    all(iy >= 1L), all(iy < n_grid)
  )

  source_weight <- rep(1 / length(source_x), length(source_x))
  mass <- numeric(n_grid * n_grid)
  mass <- add_mass(
    mass,
    ix + (iy - 1L) * n_grid,
    source_weight * (1 - fraction_x) * (1 - fraction_y)
  )
  mass <- add_mass(
    mass,
    ix + 1L + (iy - 1L) * n_grid,
    source_weight * fraction_x * (1 - fraction_y)
  )
  mass <- add_mass(
    mass,
    ix + iy * n_grid,
    source_weight * (1 - fraction_x) * fraction_y
  )
  mass <- add_mass(
    mass,
    ix + 1L + iy * n_grid,
    source_weight * fraction_x * fraction_y
  )

  matrix(mass, nrow = n_grid, ncol = n_grid)
}

# Convolve the gridded masses with the N(0, sigma_q^2 I_2) density. Zero
# padding prevents circular wraparound. Since mass contains probability masses
# and the kernel contains density values, no grid-spacing factor is introduced.
convolve_with_gaussian <- function(mass, spacing, standard_deviation) {
  n_grid <- nrow(mass)
  padded_size <- 2^ceiling(log2(2 * n_grid - 1))

  padded_mass <- matrix(0, nrow = padded_size, ncol = padded_size)
  padded_mass[seq_len(n_grid), seq_len(n_grid)] <- mass

  index <- 0:(padded_size - 1L)
  signed_offset <- ifelse(
    index <= padded_size / 2,
    index,
    index - padded_size
  ) * spacing
  gaussian_kernel <- outer(
    signed_offset,
    signed_offset,
    function(dx, dy) {
      exp(-(dx^2 + dy^2) / (2 * standard_deviation^2)) /
        (2 * pi * standard_deviation^2)
    }
  )

  convolution <- Re(
    fft(fft(padded_mass) * fft(gaussian_kernel), inverse = TRUE)
  ) / padded_size^2
  pmax(convolution[seq_len(n_grid), seq_len(n_grid)], 0)
}

compute_alternative_density <- function(q, grid, quadrature_size) {
  parameters <- curve_parameters(q)
  spacing <- grid[2L] - grid[1L]
  mass <- deposit_curve_distribution(parameters, grid, quadrature_size)
  source_mass <- sum(mass)
  density <- convolve_with_gaussian(mass, spacing, parameters$sigma_q)

  list(
    density = density,
    parameters = parameters,
    source_mass = source_mass
  )
}

# Stable evaluation of the benchmark noisy-circle density. The scaled Bessel
# function absorbs exp(-r / sigma_circle^2), leaving a well-conditioned factor.
compute_target_density <- function(grid) {
  radius <- sqrt(outer(grid^2, grid^2, "+"))
  exp(-(radius - 1)^2 / (2 * sigma_circle^2)) *
    besselI(
      radius / sigma_circle^2,
      nu = 0,
      expon.scaled = TRUE
    ) /
    (2 * pi * sigma_circle^2)
}

trapezoid_mass <- function(density, spacing) {
  edge_weight <- rep(1, nrow(density))
  edge_weight[c(1L, length(edge_weight))] <- 0.5
  sum(density * outer(edge_weight, edge_weight)) * spacing^2
}

grid <- seq(-plot_limit, plot_limit, length.out = grid_size)
grid_spacing <- grid[2L] - grid[1L]
parameters <- lapply(seq_len(9L), curve_parameters)
density_results <- vector("list", 9L)

for (q in seq_len(9L)) {
  message(sprintf("Computing p_X_%d ...", q))
  density_results[[q]] <- compute_alternative_density(
    q,
    grid,
    radial_quadrature_size
  )
}
target_density <- compute_target_density(grid)

# Numerical invariants catch missing mixture factors, FFT scaling errors, and
# accidental loss of mass at the grid-deposition stage.
source_masses <- vapply(density_results, `[[`, numeric(1), "source_mass")
grid_masses <- vapply(
  density_results,
  function(result) trapezoid_mass(result$density, grid_spacing),
  numeric(1)
)
target_grid_mass <- trapezoid_mass(target_density, grid_spacing)
stopifnot(
  max(abs(source_masses - 1)) < 1e-12,
  min(grid_masses) > 0.999,
  max(grid_masses) < 1.001,
  target_grid_mass > 0.999,
  target_grid_mass < 1.001,
  max(vapply(density_results, function(result) max(result$density), numeric(1))) <
    density_ceiling,
  max(target_density) < density_ceiling
)

validation <- data.frame(
  distribution = c(sprintf("X_%d", seq_len(9L)), "X_star"),
  integrated_mass = c(grid_masses, target_grid_mass),
  maximum_density = c(
    vapply(density_results, function(result) max(result$density), numeric(1)),
    max(target_density)
  )
)
print(validation, row.names = FALSE, digits = 8)

density_palette <- hcl.colors(256L, palette = "Inferno")
density_breaks <- seq(0, density_ceiling, length.out = 257L)
axis_ticks <- c(-1, 0, 1)
unit_theta <- seq(0, 2 * pi, length.out = 9001L)
unit_x <- cos(unit_theta)
unit_y <- sin(unit_theta)

draw_axes <- function() {
  axis(
    1,
    at = axis_ticks,
    labels = axis_ticks,
    tck = -0.025,
    cex.axis = 0.85,
    mgp = c(0, 0.22, 0)
  )
  axis(
    2,
    at = axis_ticks,
    labels = axis_ticks,
    las = 1,
    tck = -0.025,
    cex.axis = 0.85,
    mgp = c(0, 0.22, 0)
  )
  box(lwd = 0.55, col = "grey35")
}

density_title <- function(q = NULL) {
  if (is.null(q)) {
    if (isTRUE(l10n_info()[["UTF-8"]])) {
      expression(p[bolditalic(X)["\u22c6"]](bolditalic(x)))
    } else {
      expression(p[bolditalic(X)["*"]](bolditalic(x)))
    }
  } else {
    bquote(plain("PDF of")~bolditalic(X)[.(q)])
  }
}

manifold_title <- function(q = NULL) {
  if (isTRUE(l10n_info()[["UTF-8"]])) {
    if (is.null(q)) {
      expression("\u2133"["\u22c6"])
    } else {
      bquote("\u2133"[.(q)])
    }
  } else if (is.null(q)) {
    expression(italic(M)["*"])
  } else {
    bquote(italic(M)[.(q)])
  }
}

draw_panel_heading <- function(main, panel_label) {
  title(main = main, line = 1.05, cex.main = 1.35)
  mtext(
    panel_label,
    side = 3,
    line = 1.05,
    adj = 0,
    cex = 1.00,
    font = 2
  )
}

draw_density_panel <- function(density, q = NULL, panel_label) {
  image(
    grid,
    grid,
    density,
    xlim = c(-plot_limit, plot_limit),
    ylim = c(-plot_limit, plot_limit),
    zlim = c(0, density_ceiling),
    breaks = density_breaks,
    col = density_palette,
    useRaster = TRUE,
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  contour(
    grid,
    grid,
    density,
    levels = c(0.1, 0.3, 0.5),
    add = TRUE,
    drawlabels = FALSE,
    col = grDevices::adjustcolor("white", alpha.f = 0.55),
    lwd = 0.38
  )
  draw_axes()
  draw_panel_heading(density_title(q), panel_label)
}

draw_manifold_panel <- function(q = NULL, panel_label) {
  plot.new()
  plot.window(
    xlim = c(-plot_limit, plot_limit),
    ylim = c(-plot_limit, plot_limit),
    asp = 1,
    xaxs = "i",
    yaxs = "i"
  )
  abline(
    h = 0,
    v = 0,
    col = grDevices::adjustcolor("grey70", alpha.f = 0.55),
    lwd = 0.35
  )

  if (is.null(q)) {
    lines(unit_x, unit_y, col = "#145A86", lwd = 1.05)
  } else {
    p <- parameters[[q]]
    theta <- seq(
      0,
      2 * pi,
      length.out = max(9001L, 100L * p$kappa + 1L)
    )
    radius <- p$lower_radius +
      p$radial_width * (1 - cos(p$kappa * theta)) / 2
    lines(
      unit_x,
      unit_y,
      col = "grey68",
      lty = 3,
      lwd = 0.55
    )
    lines(
      radius * cos(theta),
      radius * sin(theta),
      col = "#145A86",
      lwd = 0.58
    )
  }

  draw_axes()
  draw_panel_heading(manifold_title(q), panel_label)
}

bitmap_type <- if (identical(Sys.info()[["sysname"]], "Darwin")) {
  "quartz"
} else if (capabilities("cairo")) {
  "cairo"
} else {
  getOption("bitmapType")
}

figure_family <- "sans"
if (identical(Sys.info()[["sysname"]], "Darwin")) {
  quartzFonts(
    FigureMath = quartzFont(
      c(
        "STIXTwoMath-Regular",
        "STIXGeneral-Bold",
        "STIXGeneral-Italic",
        "STIXGeneral-BoldItalic"
      )
    )
  )
  figure_family <- "FigureMath"
}

png(
  filename = output_file,
  width = 6000,
  height = 4800,
  res = 400,
  bg = "white",
  type = bitmap_type
)
par(
  mfrow = c(4, 5),
  mar = c(1.35, 1.60, 3.00, 0.35),
  oma = c(0.25, 0.25, 0.50, 0.25),
  xaxs = "i",
  yaxs = "i",
  pty = "s",
  family = figure_family
)

# Panels are lettered in row-major order. Each density is directly above its
# corresponding latent-support manifold.
panel_labels <- sprintf("(%s)", letters[seq_len(20L)])
for (q in 1:5) {
  draw_density_panel(density_results[[q]]$density, q, panel_labels[q])
}
for (q in 1:5) {
  draw_manifold_panel(q, panel_labels[5L + q])
}
for (q in 6:9) {
  draw_density_panel(
    density_results[[q]]$density,
    q,
    panel_labels[5L + q]
  )
}
draw_density_panel(target_density, panel_label = panel_labels[15L])
for (q in 6:9) {
  draw_manifold_panel(q, panel_labels[10L + q])
}
draw_manifold_panel(panel_label = panel_labels[20L])

invisible(dev.off())
message(sprintf("Wrote %s", normalizePath(output_file)))
message(
  "Integrated masses: ",
  paste(format(c(grid_masses, target_grid_mass), digits = 7), collapse = ", ")
)
message("Common density scale: [0, ", density_ceiling, "]")
