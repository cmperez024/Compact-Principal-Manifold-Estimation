## Run the fixed-design experiment with response-error SD 0.04.
find_project_root <- function(start = getwd()) {
  here <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(here, "CompactPME", "R")) &&
        dir.exists(file.path(here, "simulation_reproduction"))) return(here)
    parent <- dirname(here)
    if (identical(parent, here)) stop("Could not locate project root.")
    here <- parent
  }
}
project_root <- find_project_root()
setwd(project_root)
experiment_root <- file.path(
  project_root, "simulation_reproduction", "01_main_text",
  "04_downstream_tuning", "code_and_results"
)
Sys.setenv(
  PME_FIXED_DESIGN_EXPERIMENT =
    "1d_sd04",
  PME_RESPONSE_NOISE_SD = "0.04"
)
source(file.path(
  experiment_root, "1d_base",
  "generate-results.R"
))
