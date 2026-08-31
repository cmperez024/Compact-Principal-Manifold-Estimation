## Data-only comparison of sampling fixes for the two-dimensional flower.
## No PME or regression fit is computed.

options(stringsAsFactors = FALSE)
set.seed(20260808)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Supply the output HTML-fragment path.")
html_path <- normalizePath(args[1], winslash = "/", mustWork = FALSE)

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
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

local_library <- file.path(repo_root, "packages")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
cart2sph <- pracma::cart2sph
sph2cart <- pracma::sph2cart
source(file.path(manuscript_section, "shared", "datasets.R"))

original_phase_map <- function(z) {
  q <- (0.5 + 0.45 * sin(
    pi * z[, 1] + 0.35 * sin(pi * z[, 2])
  )) %% 1
  vertical <- 0.85 * sin(
    pi * z[, 2] + 0.25 * sin(pi * z[, 1])
  )
  angle <- 2 * pi * q
  cbind(
    sqrt(1 - vertical^2) * cos(angle),
    sqrt(1 - vertical^2) * sin(angle),
    vertical
  )
}

## If Z is uniform on [-1,1]^2, q is uniform conditionally on u2 because
## 2*u1 modulo one is uniform. The triangular-wave vertical coordinate is
## uniform on [-1,1]. The map is many-to-one in both coordinates, so Z is not
## the spherical projection index even though the induced template sample is
## uniform on S^2.
uniform_noninjective_phase_map <- function(z) {
  u1 <- (z[, 1] + 1) / 2
  u2 <- (z[, 2] + 1) / 2
  q <- (2 * u1 + 0.22 * sin(2 * pi * u2)) %% 1
  vertical <- (2 / pi) * asin(sin(2 * pi * u2))
  angle <- 2 * pi * q
  cbind(
    sqrt(1 - vertical^2) * cos(angle),
    sqrt(1 - vertical^2) * sin(angle),
    vertical
  )
}

make_data <- function(z, phase_map, noise) {
  template <- phase_map(z)
  clean <- flower2d3D_func(
    template, r0 = 1, a = 0.30, petals = 5, b = 0.50
  )
  list(
    z = z,
    template = template,
    clean = clean,
    observed = clean + noise
  )
}

n_small <- 400L
n_large <- 2000L
noise_sd <- 0.04

z_original_large <- cbind(
  z1 = runif(n_large, -1, 1),
  z2 = runif(n_large, -1, 1)
)
noise_original_large <- matrix(
  rnorm(3L * n_large, sd = noise_sd), nrow = n_large, ncol = 3L
)
z_uniform_large <- cbind(
  z1 = runif(n_large, -1, 1),
  z2 = runif(n_large, -1, 1)
)
noise_uniform_large <- matrix(
  rnorm(3L * n_large, sd = noise_sd), nrow = n_large, ncol = 3L
)

scenarios <- list(
  baseline_400 = make_data(
    z_original_large[1:n_small, , drop = FALSE],
    original_phase_map,
    noise_original_large[1:n_small, , drop = FALSE]
  ),
  original_2000 = make_data(
    z_original_large,
    original_phase_map,
    noise_original_large
  ),
  uniform_400 = make_data(
    z_uniform_large[1:n_small, , drop = FALSE],
    uniform_noninjective_phase_map,
    noise_uniform_large[1:n_small, , drop = FALSE]
  ),
  uniform_2000 = make_data(
    z_uniform_large,
    uniform_noninjective_phase_map,
    noise_uniform_large
  )
)

scenario_labels <- c(
  baseline_400 = "Original phase map, N = 400",
  original_2000 = "Original phase map, N = 2000",
  uniform_400 = "Uniform non-injective phase map, N = 400",
  uniform_2000 = "Uniform non-injective phase map, N = 2000"
)

coverage_summary <- do.call(rbind, lapply(names(scenarios), function(name) {
  template <- scenarios[[name]]$template
  longitude <- atan2(template[, 2], template[, 1]) %% (2 * pi)
  vertical <- template[, 3]
  longitude_bin <- cut(
    longitude, breaks = seq(0, 2 * pi, length.out = 13L),
    include.lowest = TRUE
  )
  vertical_bin <- cut(
    vertical, breaks = seq(-1, 1, length.out = 7L),
    include.lowest = TRUE
  )
  count <- as.numeric(table(longitude_bin, vertical_bin))
  data.frame(
    scenario = name,
    label = scenario_labels[name],
    n = nrow(template),
    equal_area_bin_cv = sd(count) / mean(count),
    empty_equal_area_bins = sum(count == 0),
    min_vertical = min(vertical),
    max_vertical = max(vertical)
  )
}))
rownames(coverage_summary) <- NULL

write.csv(
  coverage_summary,
  file.path(result_dir, "sampling-fixes-coverage-summary.csv"),
  row.names = FALSE
)
saveRDS(list(
  specification = list(
    n_small = n_small,
    n_large = n_large,
    noise_sd = noise_sd,
    original_phase_map = paste(
      "q=(0.5+0.45*sin(pi*z1+0.35*sin(pi*z2))) mod 1;",
      "s=0.85*sin(pi*z2+0.25*sin(pi*z1))"
    ),
    uniform_phase_map = paste(
      "q=(2*u1+0.22*sin(2*pi*u2)) mod 1;",
      "s=(2/pi)*asin(sin(2*pi*u2))"
    )
  ),
  scenarios = scenarios,
  coverage_summary = coverage_summary
), file.path(result_dir, "sampling-fixes-data-objects.rds"))

camera_rotation <- function(azimuth = 35, elevation = 23) {
  az <- azimuth * pi / 180
  el <- elevation * pi / 180
  rz <- matrix(c(cos(az), -sin(az), 0,
                 sin(az),  cos(az), 0,
                 0, 0, 1), 3, 3, byrow = TRUE)
  rx <- matrix(c(1, 0, 0,
                 0, cos(el), -sin(el),
                 0, sin(el),  cos(el)), 3, 3, byrow = TRUE)
  rx %*% rz
}

draw_cloud <- function(points, main, color, common_limit, point_cex) {
  rotation <- camera_rotation()
  rotated <- t(rotation %*% t(points))
  point_order <- order(rotated[, 3])
  plot(NA, xlim = c(-common_limit, common_limit),
       ylim = c(-common_limit, common_limit), asp = 1,
       axes = FALSE, xlab = "", ylab = "", main = main)
  points(rotated[point_order, 1], rotated[point_order, 2],
         pch = 16, cex = point_cex,
         col = adjustcolor(color, 0.58))
  box()
}

all_observed <- do.call(rbind, lapply(scenarios, `[[`, "observed"))
all_rotated <- t(camera_rotation() %*% t(all_observed))
flower_limit <- max(abs(all_rotated[, 1:2])) * 1.03

png(file.path(result_dir, "07-data-sampling-fixes-flower.png"),
    width = 2600, height = 2400, res = 210)
par(mfrow = c(2, 2), mar = c(1.0, 1.0, 3.8, 1.0))
colors <- c("#D55E00", "#D55E00", "#0072B2", "#0072B2")
for (i in seq_along(scenarios)) {
  name <- names(scenarios)[i]
  draw_cloud(
    scenarios[[i]]$observed,
    scenario_labels[name], colors[i], flower_limit,
    if (nrow(scenarios[[i]]$observed) <= n_small) 0.62 else 0.28
  )
}
dev.off()

template_limit <- 1.04
png(file.path(result_dir, "08-data-sampling-fixes-template.png"),
    width = 2600, height = 2400, res = 210)
par(mfrow = c(2, 2), mar = c(1.0, 1.0, 3.8, 1.0))
for (i in seq_along(scenarios)) {
  name <- names(scenarios)[i]
  draw_cloud(
    scenarios[[i]]$template,
    paste0(scenario_labels[name], "\nlatent S2 locations"),
    colors[i], template_limit,
    if (nrow(scenarios[[i]]$template) <= n_small) 0.62 else 0.28
  )
}
dev.off()

## Interactive two-scene comparison -----------------------------------------
as_js_vector <- function(x, digits = 6L) {
  paste0("[", paste(formatC(x, digits = digits, format = "fg"),
                    collapse = ","), "]")
}

js_string_vector <- function(x) {
  escaped <- gsub("\\\\", "\\\\\\\\", x)
  escaped <- gsub('"', '\\\\"', escaped)
  paste0("[", paste0('"', escaped, '"', collapse = ","), "]")
}

scenario_js <- vapply(names(scenarios), function(name) {
  value <- scenarios[[name]]
  hover <- sprintf(
    "%s<br>ID %d<br>Z = (%.3f, %.3f)",
    scenario_labels[name], seq_len(nrow(value$z)),
    value$z[, 1], value$z[, 2]
  )
  paste0(
    "{label:", js_string_vector(scenario_labels[name]), "[0]",
    ",flower:{x:", as_js_vector(value$observed[, 1]),
    ",y:", as_js_vector(value$observed[, 2]),
    ",z:", as_js_vector(value$observed[, 3]),
    ",text:", js_string_vector(hover), "}",
    ",template:{x:", as_js_vector(value$template[, 1]),
    ",y:", as_js_vector(value$template[, 2]),
    ",z:", as_js_vector(value$template[, 3]),
    ",text:", js_string_vector(hover), "}}"
  )
}, character(1))

coverage_text <- sprintf(
  "%s: bin CV %.2f; empty bins %d",
  coverage_summary$label,
  coverage_summary$equal_area_bin_cv,
  coverage_summary$empty_equal_area_bins
)

lines <- c(
  '<script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>',
  '<div id="sampling-fixes-root">',
  '  <h2>Sampling fixes for the two-dimensional flower</h2>',
  '  <div class="viz-controls" aria-label="Sampling scenario">',
  '    <button type="button" class="btn btn-primary" data-scenario="0">Original N = 400</button>',
  '    <button type="button" class="btn" data-scenario="1">Original N = 2000</button>',
  '    <button type="button" class="btn" data-scenario="2">Uniform N = 400</button>',
  '    <button type="button" class="btn" data-scenario="3">Uniform N = 2000</button>',
  '  </div>',
  '  <div id="sampling-fixes-status" class="text-small text-muted"></div>',
  '  <div id="sampling-fixes-plot" role="img" aria-label="Interactive comparison of observed two-dimensional flower point clouds and their latent spherical sampling locations."></div>',
  '</div>',
  '<style>',
  '#sampling-fixes-root { width: 100%; color: var(--foreground); }',
  '#sampling-fixes-root h2 { margin: 0 0 0.35rem 0; font-weight: 500; }',
  '#sampling-fixes-root .viz-controls { margin-bottom: 0.35rem; }',
  '#sampling-fixes-status { min-height: 1.25rem; }',
  '#sampling-fixes-plot { width: 100%; min-height: 590px; }',
  '@media (max-width: 560px) { #sampling-fixes-plot { min-height: 760px; } }',
  '</style>',
  '<script>',
  '(() => {',
  '  const root = document.getElementById("sampling-fixes-root");',
  '  const plot = document.getElementById("sampling-fixes-plot");',
  '  const status = document.getElementById("sampling-fixes-status");',
  paste0('  const scenarios = [', paste(scenario_js, collapse = ","), '];'),
  paste0('  const coverage = ', js_string_vector(coverage_text), ';'),
  '  const cssColor = name => getComputedStyle(root).getPropertyValue(name).trim();',
  '  let active = 0;',
  '  const draw = () => {',
  '    const fg = cssColor("--foreground");',
  '    const grid = cssColor("--border");',
  '    const muted = cssColor("--muted-foreground");',
  '    const flowerColor = cssColor("--viz-series-1");',
  '    const templateColor = cssColor("--viz-series-2");',
  '    const current = scenarios[active];',
  '    const size = current.flower.x.length <= 400 ? 3.2 : 1.9;',
  '    const axis = title => ({title:{text:title,font:{color:fg}},tickfont:{color:muted},gridcolor:grid,zerolinecolor:grid,showbackground:false});',
  '    const traces = [',
  '      {type:"scatter3d",mode:"markers",scene:"scene",x:current.flower.x,y:current.flower.y,z:current.flower.z,text:current.flower.text,hovertemplate:"%{text}<extra></extra>",name:"Observed flower data",marker:{size:size,color:flowerColor,opacity:0.68}},',
  '      {type:"scatter3d",mode:"markers",scene:"scene2",x:current.template.x,y:current.template.y,z:current.template.z,text:current.template.text,hovertemplate:"%{text}<extra></extra>",name:"Latent S2 locations",marker:{size:size,color:templateColor,opacity:0.68}}',
  '    ];',
  '    const narrow = root.getBoundingClientRect().width < 560;',
  '    const scene1Domain = narrow ? {x:[0,1],y:[0.52,1]} : {x:[0,0.49],y:[0,1]};',
  '    const scene2Domain = narrow ? {x:[0,1],y:[0,0.48]} : {x:[0.51,1],y:[0,1]};',
  '    const layout = {margin:{l:0,r:0,t:36,b:12},paper_bgcolor:"rgba(0,0,0,0)",plot_bgcolor:"rgba(0,0,0,0)",font:{color:fg},showlegend:false,',
  '      annotations:[{text:"Observed flower data",x:narrow?0.5:0.245,y:narrow?1:1.03,xref:"paper",yref:"paper",showarrow:false,font:{color:fg,size:14}},{text:"Latent S2 locations",x:narrow?0.5:0.755,y:narrow?0.48:1.03,xref:"paper",yref:"paper",showarrow:false,font:{color:fg,size:14}}],',
  '      scene:{domain:scene1Domain,aspectmode:"data",xaxis:axis("X1"),yaxis:axis("X2"),zaxis:axis("X3"),camera:{eye:{x:1.45,y:1.45,z:0.95}}},',
  '      scene2:{domain:scene2Domain,aspectmode:"cube",xaxis:axis("T1"),yaxis:axis("T2"),zaxis:axis("T3"),camera:{eye:{x:1.45,y:1.45,z:0.95}}}};',
  '    status.textContent = coverage[active];',
  '    Plotly.react(plot,traces,layout,{responsive:true,displaylogo:false,scrollZoom:true});',
  '  };',
  '  root.querySelectorAll("button[data-scenario]").forEach(button => {',
  '    button.addEventListener("click", () => {',
  '      active = Number(button.dataset.scenario);',
  '      root.querySelectorAll("button[data-scenario]").forEach((peer,index) => {',
  '        peer.classList.toggle("btn-primary", index === active);',
  '      });',
  '      draw();',
  '    });',
  '  });',
  '  draw();',
  '  new ResizeObserver(draw).observe(root);',
  '  new MutationObserver(draw).observe(document.documentElement,{attributes:true,attributeFilter:["class","style"]});',
  '})();',
  '</script>'
)

dir.create(dirname(html_path), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, html_path, useBytes = TRUE)

print(coverage_summary)
cat(html_path, "\n")
