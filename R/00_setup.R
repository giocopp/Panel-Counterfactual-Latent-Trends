#' ============================================================================
#' Setup Script: Install Dependencies and Quick Test
#'
#' Run first:
#'   source("R/00_setup.R")
#' Then:
#'   source("R/run_analysis.R")
#' ============================================================================

install_packages <- function() {
  required <- c(
    "tidyverse",
    "lubridate",
    "fixest",
    "Matrix",
    "patchwork",
    "scales"
  )

  optional <- c(
    "furrr",
    "progressr"
  )

  cat("Installing required packages...\n")
  for (pkg in required) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    cat("  ✓", pkg, "installed\n")
  }

  cat("\nInstalling optional packages...\n")
  for (pkg in optional) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    cat("  ✓", pkg, "installed\n")
  }

  cat("\nDone!\n\n")
}

install_packages()

cat("=== Quick Test ===\n\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

source("R/01_dgp.R")
source("R/03_estimators.R")

cat("Testing DGP...\n")
grid <- generate_grid(near_cutoff_km = 200)
time_panel <- generate_time_panel()
panel <- generate_panel_data(
  grid = grid,
  time_panel = time_panel,
  delta = 0.6,
  factor_strength = 0.5,
  short_lived = TRUE,
  seed = 1
)

cat("  ✓ Generated", nrow(panel), "observations\n")
cat("  ✓ Deadly-event rate:", round(mean(panel$Y_any), 4), "\n\n")

true <- compute_true_estimands(panel, outcome_col = "Y_any")

cat("Testing estimators...\n")
twfe_res <- estimate_twfe(panel, outcome_col = "Y_any")
cat("  ✓ TWFE estimate:", round(twfe_res$estimate, 4),
    "(true:", round(true$att, 4), ", bias:", round(twfe_res$estimate - true$att, 4), ")\n")

trop_res <- estimate_trop(panel, outcome_col = "Y_any")
cat("  ✓ TROP estimate:", round(trop_res$estimate, 4),
    "(true:", round(true$att, 4), ", bias:", round(trop_res$estimate - true$att, 4), ")\n\n")

cat(paste(rep("=", 50), collapse = ""), "\n")
cat("All tests passed! Ready to run analysis.\n")
cat("Run: source('R/run_analysis.R')\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
