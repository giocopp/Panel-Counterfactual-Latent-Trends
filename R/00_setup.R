#' ============================================================================
#' Setup Script: Install Dependencies and Quick Test
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
    cat("  ✓", pkg, "\n")
  }

  # Ensure remotes is available
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }

  # synthdid from GitHub
  cat("\nInstalling synthdid from GitHub...\n")
  if (!requireNamespace("synthdid", quietly = TRUE)) {
    tryCatch(
      {
        remotes::install_github("synth-inference/synthdid")
        cat("  ✓ synthdid\n")
      },
      error = function(e) {
        cat("  ✗ synthdid failed - will use custom implementation\n")
      }
    )
  } else {
    cat("  ✓ synthdid (already installed)\n")
  }

  # Skip MCPanel - won't compile on modern macOS
  cat("\nNote: MCPanel skipped (C++ compilation issues on macOS).\n")
  cat("      Using custom matrix completion implementation.\n")

  cat("\nInstalling optional packages...\n")
  for (pkg in optional) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    cat("  ✓", pkg, "\n")
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
cat(
  "  ✓ TWFE estimate:", round(twfe_res$estimate, 4),
  "(bias:", round(twfe_res$estimate - true$att, 4), ")\n"
)

mc_res <- estimate_mc(panel, outcome_col = "Y_any")
cat(
  "  ✓ MC estimate:", round(mc_res$estimate, 4),
  "(bias:", round(mc_res$estimate - true$att, 4), ")\n"
)

sdid_res <- estimate_sdid(panel, outcome_col = "Y_any")
cat(
  "  ✓ SDID estimate:", round(sdid_res$estimate, 4),
  "(bias:", round(sdid_res$estimate - true$att, 4), ")\n"
)

trop_res <- estimate_trop(panel, outcome_col = "Y_any")
cat(
  "  ✓ TROP estimate:", round(trop_res$estimate, 4),
  "(bias:", round(trop_res$estimate - true$att, 4), ")\n\n"
)

cat("True ATT:", round(true$att, 4), "\n\n")

cat(paste(rep("=", 50), collapse = ""), "\n")
cat("All tests passed! Ready to run analysis.\n")
cat("Run: source('R/run_analysis.R')\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
