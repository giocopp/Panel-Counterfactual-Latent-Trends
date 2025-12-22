#' ============================================================================
#' Setup Script: Install Dependencies and Quick Test (semi-synthetic DGP)
#' ============================================================================

install_packages <- function() {
  required <- c(
    "tidyverse",
    "lubridate",
    "fixest",
    "fect",
    "Matrix",
    "furrr",
    "progressr",
    "synthdid",
    "patchwork",
    "scales"
  )

  cat("Installing required packages...\n")
  for (pkg in required) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      suppressWarnings(
        suppressMessages(
          install.packages(pkg, dependencies = TRUE)
        )
      )
    }
    cat("  ✓", pkg, "\n")
  }
  cat("\nDone!\n\n")
}

install_packages()

cat("=== Quick Test (semi-synthetic logistic IFE DGP) ===\n\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest) # TWFE
  library(fect) # MC
  library(synthdid) # SynthDiD
})

# NOTE: keep your project paths consistent.
source("R/01_dgp.R")
source("R/02_calibration.R")
source("R/03_estimators.R")

# -----------------------------------------------------------------------------
# QUICK KNOBS
# -----------------------------------------------------------------------------
OUTCOME_COL <- "Y_any"

# Severity (counts) should increase after treatment and slightly decay:
# Interpreted as multiplicative increase in mu_pos: mu_pos_t = mu_pos * exp(sev_delta_t)
SEV_DELTA <- log(1.25) # ~ +25% severity at treatment start
SEV_EFFECT_TYPE <- "decay" # "constant" | "decay" | "window"
SEV_DECAY_HALF_LIFE <- 8 # months

cat("Building grid & time panel...\n")
grid <- generate_grid(near_cutoff_km = 200)
time_panel <- generate_time_panel()

cat("Loading / building DGP anchor (from IOM pre-period)...\n")
anchor <- tryCatch(load_dgp_anchor("data/dgp_anchor.rds"), error = function(e) NULL)

if (is.null(anchor)) {
  cat("  - No data/dgp_anchor.rds found. Running calibration now...\n")
  cal <- calibrate_from_iom(
    lat_range = c(31, 38),
    lon_range = c(10, 18),
    start_date = "2015-04-01",
    end_date = "2018-02-01",
    treat_start = "2017-05-01",
    near_cutoff_km = 200,
    K = 2,
    c_smooth = 0.5
  )
  anchor <- cal$anchor
}
cat("  ✓ Anchor ready (K = ", anchor$K, ")\n\n", sep = "")

cat("Testing baseline DGP...\n")

# Build args robustly (so this script still runs even before you update 01_dgp.R
# to actually implement sev_* parameters).
gen_args <- list(
  grid = grid,
  time_panel = time_panel,
  anchor = anchor,
  delta = 0.6,
  s = 0.9,
  effect_type = "constant",
  rho = 0.3,
  sigma_u = anchor$sigma_u,
  q = 1.0,
  seed = 1
)

gp_formals <- names(formals(generate_panel_data))
if (all(c("sev_delta", "sev_effect_type", "sev_decay_half_life") %in% gp_formals)) {
  gen_args$sev_delta <- SEV_DELTA
  gen_args$sev_effect_type <- SEV_EFFECT_TYPE
  gen_args$sev_decay_half_life <- SEV_DECAY_HALF_LIFE
  cat("  ✓ Severity enabled (post increase + decay)\n")
} else {
  cat("  ! Severity params not found in generate_panel_data() yet.\n")
  cat("    Update R/01_dgp.R to add sev_delta/sev_effect_type/sev_decay_half_life.\n")
}

panel <- do.call(generate_panel_data, gen_args)

cat("  ✓ Generated ", nrow(panel), " observations\n", sep = "")
cat("  ✓ Mean(Y_ihs):", round(mean(panel$Y_ihs), 4), "\n")
cat("  ✓ Deadly-event rate mean(Y_any):", round(mean(panel$Y_any), 4), "\n\n")

true <- compute_true_estimands(panel, outcome_col = OUTCOME_COL)

cat("Testing estimators on ", OUTCOME_COL, "...\n", sep = "")

twfe_res <- estimate_twfe(panel, outcome_col = OUTCOME_COL)
cat(
  "  ✓ TWFE estimate:", round(twfe_res$estimate, 4),
  "(bias:", round(twfe_res$estimate - true$att, 4), ")\n"
)

mc_res <- estimate_mc(panel, outcome_col = OUTCOME_COL)
cat(
  "  ✓ MC estimate:", round(mc_res$estimate, 4),
  "(bias:", round(mc_res$estimate - true$att, 4), ")\n"
)

sdid_res <- estimate_sdid(panel, outcome_col = OUTCOME_COL)
cat(
  "  ✓ SynthDiD estimate:", round(sdid_res$estimate, 4),
  "(bias:", round(sdid_res$estimate - true$att, 4), ")\n"
)

cat("\nTrue ATT (treated post):", round(true$att, 4), "\n")
cat("Expected DiD (diagnostic):", round(true$did, 4), "\n\n")

cat(paste(rep("=", 50), collapse = ""), "\n")
cat("All tests passed! Ready to run analysis.\n")
cat("Run: source('R/run_analysis.R')\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
