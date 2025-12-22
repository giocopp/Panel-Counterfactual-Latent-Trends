#' ============================================================================
#' Main Analysis Script (End-to-End)
#'
#' Theme: Panel counterfactual estimators for policy shocks under latent trends
#' Application-motivated DGP: Italy–Libya MoU (May 2017 operational onset)
#'
#' Estimators compared:
#'   - TWFE: Two-way fixed effects (fixest)
#'   - Matrix Completion: Low-rank imputation 
#'   - TROP: Triply robust panel estimator (custom)
#'   - SynthDiD: Synthetic Difference-in-Differences (synthdid)
#'
#' Outputs used by writeup.qmd:
#'   - output/power_summary.csv
#'   - output/scenario_summary.csv
#'   - figures/power_curves.pdf
#'   - figures/bias_curve.pdf
#'   - figures/coverage_curve.pdf
#' ============================================================================

# Set seed for reproducibility
set.seed(42)

# -----------------------------------------------------------------------------
# Quiet install helper (keeps script self-contained + consistent with setup)
# -----------------------------------------------------------------------------
quiet_install <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      suppressWarnings(suppressMessages(
        install.packages(pkg, dependencies = TRUE)
      ))
    }
  }
}

# Packages used in this script (quietly install if missing)
quiet_install(c(
  "tidyverse", "fixest", "patchwork", "scales",
  "furrr", "progressr", "synthdid"
))

# -----------------------------------------------------------------------------
# Parallel setup (optional)
# -----------------------------------------------------------------------------
if (requireNamespace("furrr", quietly = TRUE)) {
  suppressPackageStartupMessages(library(furrr))
  n_workers <- max(1, parallel::detectCores() - 1)
  plan(multisession, workers = n_workers)
  cat("Parallel processing enabled:", n_workers, "workers\n")
}

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest) # TWFE
  library(patchwork)
  library(scales)
  library(synthdid) # SynthDiD
})

# -----------------------------------------------------------------------------
# Source project code (functions exist only after this)
# -----------------------------------------------------------------------------
source("R/01_dgp.R")
source("R/02_calibration.R")
source("R/03_estimators.R")
source("R/04_simulation.R")

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

OUTCOME_COL <- "Y_any"
START_DATE <- "2015-04-01"
END_DATE <- "2018-02-01"
TREAT_START <- "2017-05-01"
NEAR_CUTOFF <- 200

# -----------------------------------------------------------------------------
# Estimator set (treat sdid like the others: include it => must be installed)
# -----------------------------------------------------------------------------
ESTIMATORS <- c("twfe", "mc", "trop", "sdid")
cat("Estimators:", paste(ESTIMATORS, collapse = ", "), "\n\n")

# Enforce SynthDiD availability if requested
if ("sdid" %in% ESTIMATORS && !requireNamespace("synthdid", quietly = TRUE)) {
  stop("SynthDiD requires the 'synthdid' package. Install it with: install.packages('synthdid')")
}

# =============================================================================
# PART 0: LOAD CALIBRATION FROM SAVED FILE
# =============================================================================
calibration <- NULL
calibration_path <- "data/calibration_targets.rds"

if (file.exists(calibration_path)) {
  calibration <- readRDS(calibration_path)
  cat("\nCalibration loaded from:", calibration_path, "\n")
  cat("  Target pre-period event rate (Y_any):", round(calibration$target_event_rate, 4), "\n")
  cat(
    "  Calibrated baseline_mortality:", signif(calibration$baseline_mortality, 4),
    "(", calibration$calibration_method, ")\n\n"
  )
} else {
  cat("\nNo calibration file found at", calibration_path, "\n")
  cat("Using default DGP baselines.\n")
  cat("Run calibrate_from_iom() first to generate calibration.\n\n")
}

cat("\n============================================================\n")
cat("PART 1: DGP (Zambiasi-style spatial panel)\n")
cat("============================================================\n\n")

grid <- generate_grid(near_cutoff_km = NEAR_CUTOFF)

time_panel <- generate_time_panel(
  start_date = START_DATE,
  end_date = END_DATE,
  treat_start = TREAT_START,
  seasonality_index = if (!is.null(calibration)) calibration$seasonality_index else NULL
)

# Writeup-friendly calibration diagnostics (function should be safe if IOM file missing)
if (!is.null(calibration)) {
  make_calibration_figures(
    calibration = calibration,
    grid = grid,
    time_panel = time_panel,
    treat_start = TREAT_START,
    fig_dir = "figures"
  )
}

example_data <- generate_panel_data(
  grid = grid,
  time_panel = time_panel,
  baseline_mortality = if (!is.null(calibration)) calibration$baseline_mortality else 0.02,
  delta = 1.2,
  factor_strength = 0.5,
  loading_correlation = 0.5, # Correlated loadings -> violates parallel trends
  underreporting_rate = 0.0,
  short_lived = FALSE,
  seed = 123
)

true <- compute_true_estimands(example_data, outcome_col = OUTCOME_COL)

cat("Units:", n_distinct(example_data$unit_id), "\n")
cat("Months:", n_distinct(example_data$time_id), "\n")
cat("Observations:", nrow(example_data), "\n\n")

cat("Share near Libya (units):", round(mean(grid$near_libya), 3), "\n")
cat("Deadly-event rate (overall):", percent(mean(example_data$Y_any)), "\n")
cat("Deadly-event rate (pre):", percent(mean(example_data$Y_any[example_data$post_treat == 0])), "\n")
cat("Deadly-event rate (post):", percent(mean(example_data$Y_any[example_data$post_treat == 1])), "\n\n")

cat("True ATT (implied):", round(true$att, 4), "\n")
cat("Expected DiD (implied):", round(true$did, 4), "\n\n")

cat("\n============================================================\n")
cat("PART 2: One-dataset estimation example\n")
cat("============================================================\n\n")

example_results <- run_all_estimators(
  example_data,
  outcome_col = OUTCOME_COL,
  true_att = true$att,
  estimators = ESTIMATORS
)

print(example_results %>% select(method, estimate, se, true_att, bias, p_value))
write_csv(example_results, "output/example_estimates.csv")

cat("\n============================================================\n")
cat("PART 3: Power analysis\n")
cat("============================================================\n\n")

power_results <- run_power_analysis(
  baseline_mortality = if (!is.null(calibration)) calibration$baseline_mortality else 0.02,
  seasonality_index = if (!is.null(calibration)) calibration$seasonality_index else NULL,
  delta_values = c(0, 0.4, 0.8, 1.2, 1.6), 
  n_sims = 200,
  estimators = ESTIMATORS,
  outcome_col = OUTCOME_COL,
  base_params = list(
    factor_strength = 0.5,
    loading_correlation = 0.5,
    underreporting_rate = 0.0,
    short_lived = FALSE
  )
)

power_summary <- summarize_power_results(power_results)
write_csv(power_summary, "output/power_summary.csv")

p_power <- plot_power_curves(power_summary)
ggsave("figures/power_curves.pdf", p_power, width = 7.2, height = 4.6)

p_bias <- plot_bias(power_summary)
ggsave("figures/bias_curve.pdf", p_bias, width = 7.2, height = 4.6)

p_cov <- plot_coverage(power_summary)
ggsave("figures/coverage_curve.pdf", p_cov, width = 7.2, height = 4.6)

cat("Saved figures to figures/\n")
cat("Saved power summary to output/power_summary.csv\n\n")

cat("\n============================================================\n")
cat("PART 4: Scenario analysis\n")
cat("============================================================\n\n")

scenario_results <- run_scenario_analysis(
  n_sims = 100,
  estimators = ESTIMATORS,
  outcome_col = OUTCOME_COL
)

scenario_summary <- summarize_scenario_results(scenario_results)
write_csv(scenario_summary, "output/scenario_summary.csv")

cat("Saved scenario summary to output/scenario_summary.csv\n\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("Now render the writeup: quarto render writeup.qmd\n")
cat("============================================================\n")
