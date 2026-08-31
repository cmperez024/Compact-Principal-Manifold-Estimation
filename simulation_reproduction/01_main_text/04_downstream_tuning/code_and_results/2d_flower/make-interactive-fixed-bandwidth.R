## Interactive comparison for the fixed-bandwidth local-linear result:
## h = 0.06 and lambda = 7.19685673001151e-06.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Supply the output HTML-fragment path.")
output_path <- normalizePath(args[1], winslash = "/", mustWork = FALSE)

find_repo_root <- function(start = getwd()) {
  here <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(here, "CompactPME", "R"))) return(here)
    parent <- dirname(here)
    if (identical(parent, here)) stop("Could not locate repository root.")
    here <- parent
  }
}

repo_root <- find_repo_root()
setwd(repo_root)
manuscript_section <- file.path(repo_root, "simulation_reproduction",
  "01_main_text", "04_downstream_tuning")
notebook_dir <- file.path(
  manuscript_section, "code_and_results", "2d_flower"
)
result_dir <- file.path(notebook_dir, "results")
cache_dir <- file.path(notebook_dir, "cache")
local_library <- file.path(repo_root, "packages")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))

cart2sph <- pracma::cart2sph
sph2cart <- pracma::sph2cart
source(file.path(manuscript_section, "shared", "datasets.R"))
source(file.path(repo_root, "CompactPME", "R", "qm.R"))
source(file.path(repo_root, "CompactPME", "R", "spline2d.R"))
source(file.path(repo_root, "CompactPME", "R", "project_optimize2.R"))
ones <- pracma::ones
eye <- pracma::eye

objects <- readRDS(file.path(
  result_dir, "extrinsic-regression-2d-flower-objects.rds"
))
z_external <- objects$data$z_external
response_observed <- objects$data$response_observed
fold_id <- objects$data$fold_id
lambda_grid <- objects$tuning$lambda_grid
target_lambda <- lambda_grid[7L]
target_bandwidth <- 0.06
lambda_index <- 7L

normalize_rows <- function(x) x / sqrt(rowSums(x^2))

fit_pme_s2 <- function(x, lambda, max_iterations = 15L) {
  projection <- normalize_rows(scale(x, scale = FALSE))
  fit <- spline2d(x, projection, lambda)
  projection <- project_optimize2(x, fit, projection)
  for (iteration in 2:max_iterations) {
    fit <- spline2d(x, projection, lambda)
    projection <- project_optimize2(x, fit, projection)
  }
  list(fit = fit, projection = projection, denoised = fit(projection))
}

ambient_local_linear <- function(z_train, y_train, z_new, bandwidth) {
  output <- matrix(NA_real_, nrow = nrow(z_new), ncol = ncol(y_train))
  for (j in seq_len(nrow(z_new))) {
    delta <- sweep(z_train, 2L, z_new[j, ], "-")
    weights <- exp(-0.5 * rowSums(delta^2) / bandwidth^2)
    design <- cbind(1, delta)
    normal_matrix <- crossprod(design, weights * design)
    right_hand_side <- crossprod(design, weights * y_train)
    scale_value <- max(1, max(abs(diag(normal_matrix))))
    coefficient <- if (is.finite(rcond(normal_matrix)) &&
                       rcond(normal_matrix) > 1e-10) {
      tryCatch(solve(normal_matrix, right_hand_side),
               error = function(e) NULL)
    } else {
      NULL
    }
    if (is.null(coefficient)) {
      coefficient <- tryCatch(
        solve(normal_matrix + diag(1e-8 * scale_value, 3L),
              right_hand_side),
        error = function(e) NULL
      )
    }
    if (is.null(coefficient) || any(!is.finite(coefficient))) {
      output[j, ] <- colSums(weights * y_train) / sum(weights)
    } else {
      output[j, ] <- coefficient[1, ]
    }
  }
  output
}

surface_grid <- fibonaccisphere(1800L)
project_grid <- function(points, surface) {
  distance_squared <- outer(rowSums(points^2), rowSums(surface^2), "+") -
    2 * points %*% t(surface)
  surface[max.col(-distance_squared, ties.method = "first"), , drop = FALSE]
}

## Honest OOF predictions with h fixed at 0.06.
oof_prediction <- matrix(NA_real_, nrow(response_observed), 3L)
for (fold in 1:5) {
  training_id <- which(fold_id != fold)
  validation_id <- which(fold_id == fold)
  pme <- readRDS(file.path(
    cache_dir,
    sprintf("pme-lambda-%02d-fold-%d.rds", lambda_index, fold)
  ))
  ambient <- ambient_local_linear(
    z_external[training_id, , drop = FALSE],
    pme$denoised,
    z_external[validation_id, , drop = FALSE],
    target_bandwidth
  )
  oof_prediction[validation_id, ] <- project_grid(
    ambient, pme$fit(surface_grid)
  )
}

## Full-data surface at the same lambda for geometric display only.
full_pme <- fit_pme_s2(response_observed, target_lambda)

fold_raw_rmse <- vapply(1:5, function(fold) {
  ids <- which(fold_id == fold)
  sqrt(mean(rowSums(
    (response_observed[ids, , drop = FALSE] -
       oof_prediction[ids, , drop = FALSE])^2
  )))
}, numeric(1))
fold_oracle_rmse <- vapply(1:5, function(fold) {
  ids <- which(fold_id == fold)
  sqrt(mean(rowSums(
    (objects$data$conditional_mean[ids, , drop = FALSE] -
       oof_prediction[ids, , drop = FALSE])^2
  )))
}, numeric(1))

saveRDS(list(
  lambda = target_lambda,
  bandwidth = target_bandwidth,
  z_external = z_external,
  response_observed = response_observed,
  conditional_mean = objects$data$conditional_mean,
  fold_id = fold_id,
  oof_prediction = oof_prediction,
  raw_rmse = sqrt(mean(rowSums((response_observed - oof_prediction)^2))),
  oracle_rmse = sqrt(mean(rowSums(
    (objects$data$conditional_mean - oof_prediction)^2
  ))),
  fold_raw_rmse = fold_raw_rmse,
  fold_oracle_rmse = fold_oracle_rmse,
  mean_fold_raw_rmse = mean(fold_raw_rmse),
  mean_fold_oracle_rmse = mean(fold_oracle_rmse),
  full_pme = full_pme
), file.path(result_dir, "fixed-bandwidth-local-linear-objects.rds"))

theta <- seq(0, 2 * pi, length.out = 55L)
vertical <- seq(-0.98, 0.98, length.out = 31L)
parameter_grid <- expand.grid(theta = theta, vertical = vertical)
template_grid <- cbind(
  sqrt(1 - parameter_grid$vertical^2) * cos(parameter_grid$theta),
  sqrt(1 - parameter_grid$vertical^2) * sin(parameter_grid$theta),
  parameter_grid$vertical
)
true_values <- flower2d3D_func(
  template_grid, r0 = 1, a = 0.30, petals = 5, b = 0.50
)
pme_values <- full_pme$fit(template_grid)

as_js_vector <- function(x, digits = 6L) {
  paste0("[", paste(formatC(x, digits = digits, format = "fg"),
                    collapse = ","), "]")
}

as_js_matrix <- function(x, n_row, n_col, digits = 6L) {
  value <- matrix(x, nrow = n_row, ncol = n_col)
  rows <- apply(value, 1L, as_js_vector, digits = digits)
  paste0("[", paste(rows, collapse = ","), "]")
}

js_string_vector <- function(x) {
  escaped <- gsub("\\\\", "\\\\\\\\", x)
  escaped <- gsub('"', '\\\\"', escaped)
  paste0("[", paste0('"', escaped, '"', collapse = ","), "]")
}

n_theta <- length(theta)
n_vertical <- length(vertical)
true_x <- as_js_matrix(true_values[, 1], n_theta, n_vertical)
true_y <- as_js_matrix(true_values[, 2], n_theta, n_vertical)
true_z <- as_js_matrix(true_values[, 3], n_theta, n_vertical)
pme_x <- as_js_matrix(pme_values[, 1], n_theta, n_vertical)
pme_y <- as_js_matrix(pme_values[, 2], n_theta, n_vertical)
pme_z <- as_js_matrix(pme_values[, 3], n_theta, n_vertical)

point_error <- sqrt(rowSums((response_observed - oof_prediction)^2))
observed_hover <- sprintf(
  "Observation %d<br>Z = (%.3f, %.3f)",
  seq_len(nrow(z_external)), z_external[, 1], z_external[, 2]
)
prediction_hover <- sprintf(
  "OOF prediction %d<br>Z = (%.3f, %.3f)<br>Raw error = %.3f",
  seq_len(nrow(z_external)), z_external[, 1], z_external[, 2], point_error
)

lines <- c(
  '<script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>',
  '<div id="pme-2d-fixed-root">',
  '  <h2>Two-dimensional flower at h = 0.06 and lambda = 7.20e-6</h2>',
  '  <div class="text-small text-muted">Drag to rotate; scroll or pinch to zoom. Click legend labels to show or hide layers.</div>',
  '  <div id="pme-2d-fixed-plot" role="img" aria-label="Interactive three-dimensional comparison of noisy flower responses, the true flower surface, the fitted CPME surface, and local-linear out-of-fold predictions."></div>',
  '  <div class="sr-only">The local-linear out-of-fold raw RMSE is 0.1976. The fitted CPME surface uses lambda 7.20e-6 and the regression bandwidth is 0.06.</div>',
  '</div>',
  '<style>',
  '#pme-2d-fixed-root { width: 100%; color: var(--foreground); }',
  '#pme-2d-fixed-root h2 { margin: 0 0 0.25rem 0; font-weight: 500; }',
  '#pme-2d-fixed-plot { width: 100%; min-height: 610px; }',
  '@media (max-width: 480px) { #pme-2d-fixed-plot { min-height: 470px; } }',
  '</style>',
  '<script>',
  '(() => {',
  '  const root = document.getElementById("pme-2d-fixed-root");',
  '  const plot = document.getElementById("pme-2d-fixed-plot");',
  '  const cssColor = name => getComputedStyle(root).getPropertyValue(name).trim();',
  paste0('  const observed = {x:', as_js_vector(response_observed[, 1]),
         ',y:', as_js_vector(response_observed[, 2]),
         ',z:', as_js_vector(response_observed[, 3]),
         ',text:', js_string_vector(observed_hover), '};'),
  paste0('  const predicted = {x:', as_js_vector(oof_prediction[, 1]),
         ',y:', as_js_vector(oof_prediction[, 2]),
         ',z:', as_js_vector(oof_prediction[, 3]),
         ',text:', js_string_vector(prediction_hover), '};'),
  paste0('  const trueSurface = {x:', true_x, ',y:', true_y, ',z:', true_z, '};'),
  paste0('  const pmeSurface = {x:', pme_x, ',y:', pme_y, ',z:', pme_z, '};'),
  '  const draw = () => {',
  '    const fg = cssColor("--foreground");',
  '    const grid = cssColor("--border");',
  '    const muted = cssColor("--muted-foreground");',
  '    const trueColor = cssColor("--viz-series-1");',
  '    const pmeColor = cssColor("--viz-series-2");',
  '    const observedColor = cssColor("--viz-series-3");',
  '    const predictionColor = cssColor("--viz-series-4");',
  '    const traces = [',
  '      {type:"surface", x:trueSurface.x, y:trueSurface.y, z:trueSurface.z,',
  '       name:"True flower surface", opacity:0.16, showscale:false, hoverinfo:"skip",',
  '       colorscale:[[0,trueColor],[1,trueColor]], showlegend:true},',
  '      {type:"surface", x:pmeSurface.x, y:pmeSurface.y, z:pmeSurface.z,',
  '       name:"Full-data CPME surface", opacity:0.42, showscale:false, hoverinfo:"skip",',
  '       colorscale:[[0,pmeColor],[1,pmeColor]], showlegend:true},',
  '      {type:"scatter3d", mode:"markers", x:observed.x, y:observed.y, z:observed.z,',
  '       text:observed.text, hovertemplate:"%{text}<extra></extra>",',
  '       name:"Noisy observations", marker:{size:2.4,color:observedColor,opacity:0.32}},',
  '      {type:"scatter3d", mode:"markers", x:predicted.x, y:predicted.y, z:predicted.z,',
  '       text:predicted.text, hovertemplate:"%{text}<extra></extra>",',
  '       name:"Local-linear OOF predictions", marker:{size:3.4,color:predictionColor,opacity:0.88}}',
  '    ];',
  '    const axis = title => ({title:{text:title,font:{color:fg}},tickfont:{color:muted},',
  '      gridcolor:grid,zerolinecolor:grid,showbackground:false});',
  '    const layout = {margin:{l:0,r:0,t:12,b:44},paper_bgcolor:"rgba(0,0,0,0)",',
  '      plot_bgcolor:"rgba(0,0,0,0)",font:{color:fg},legend:{orientation:"h",x:0,y:-0.04},',
  '      scene:{aspectmode:"data",xaxis:axis("X1"),yaxis:axis("X2"),zaxis:axis("X3"),',
  '        camera:{eye:{x:1.45,y:1.45,z:0.95}}}};',
  '    Plotly.react(plot,traces,layout,{responsive:true,displaylogo:false,scrollZoom:true});',
  '  };',
  '  draw();',
  '  new ResizeObserver(() => Plotly.Plots.resize(plot)).observe(root);',
  '  new MutationObserver(draw).observe(document.documentElement,{attributes:true,attributeFilter:["class","style"]});',
  '})();',
  '</script>'
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, output_path, useBytes = TRUE)
cat(output_path, "\n")
