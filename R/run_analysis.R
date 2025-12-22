#' ============================================================================
#' Main Analysis Script
#' Outcome: Y_any (binary indicator of any deaths in cell-month)
#' Estimators: TWFE, Matrix Completion, SynthDiD
#' ============================================================================

set.seed(42)

quiet_install <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      suppressWarnings(suppressMessages(install.packages(pkg, dependencies = TRUE)))
    }
  }
}

quiet_install(c("tidyverse", "fixest", "patchwork", "scales", "furrr", "progressr", "synthdid", "fect", "lubridate"))

if (requireNamespace("furrr", quietly = TRUE)) {
  suppressPackageStartupMessages(library(furrr))
  n_workers <- max(1, parallel::detectCores() - 1)
  plan(multisession, workers = n_workers)
  cat("Parallel processing enabled:", n_workers, "workers\n")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
  library(patchwork)
  library(scales)
})

source("R/01_dgp.R")
source("R/02_calibration.R")
source("R/03_estimators.R")
source("R/04_simulation.R")

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("output", showWarnings = FALSE, recursive = TRUE)
dir.create("data", showWarnings = FALSE, recursive = TRUE)

OUTCOME_COL <- "Y_any"
START_DATE <- "2015-04-01"
END_DATE <- "2018-02-01"
TREAT_START <- "2017-05-01"
NEAR_CUTOFF <- 200
LAT_RANGE <- c(31, 38)
LON_RANGE <- c(10, 18)
K_FACTORS <- 2

ESTIMATORS <- c("twfe", "mc", "sdid")
cat("Estimators:", paste(ESTIMATORS, collapse = ", "), "\n")
cat("Outcome:", OUTCOME_COL, "\n\n")

if ("sdid" %in% ESTIMATORS && !requireNamespace("synthdid", quietly = TRUE)) {
  stop("SynthDiD requires the 'synthdid' package. Install it with: install.packages('synthdid')")
}

# =============================================================================
# PART 0: Anchor
# =============================================================================
anchor_path <- "data/dgp_anchor.rds"

anchor <- if (file.exists(anchor_path)) {
  cat("Loaded anchor from:", anchor_path, "\n")
  load_dgp_anchor(anchor_path)
} else {
  cat("Anchor not found. Building from IOM (pre-period)...\n")
  cal <- calibrate_from_iom(
    lat_range = LAT_RANGE, lon_range = LON_RANGE,
    start_date = START_DATE, end_date = END_DATE,
    treat_start = TREAT_START, near_cutoff_km = NEAR_CUTOFF,
    K = K_FACTORS, c_smooth = 0.5
  )
  if (is.null(cal) || is.null(cal$anchor)) stop("Failed to build anchor.")
  cal$anchor
}

cat("Anchor K:", anchor$K, "\n")
cat("Anchor sigma_u:", round(anchor$sigma_u, 4), "\n\n")

# =============================================================================
# PART 1: One example dataset
# =============================================================================
grid <- generate_grid(lat_range = LAT_RANGE, lon_range = LON_RANGE, near_cutoff_km = NEAR_CUTOFF)
time_panel <- generate_time_panel(start_date = START_DATE, end_date = END_DATE, treat_start = TREAT_START)

example_data <- generate_panel_data(
  grid = grid,
  time_panel = time_panel,
  anchor = anchor,
  delta = 0.6,
  s = 0.9,
  effect_type = "constant",
  rho = 0.3,
  sigma_u = anchor$sigma_u,
  q = 1.0,
  mu_pos_cap = 30,
  seed = 123
)

true <- compute_true_estimands(example_data, outcome_col = OUTCOME_COL, treat_type = "post")

cat("Units:", n_distinct(example_data$unit_id), "\n")
cat("Months:", n_distinct(example_data$time_id), "\n")
cat("Obs:", nrow(example_data), "\n")
cat("Mean(Y_any):", round(mean(example_data$Y_any), 4), "\n")
cat("True ATT:", round(true$att, 4), "\n\n")

example_results <- run_all_estimators(
  example_data,
  outcome_col = OUTCOME_COL,
  true_att = true$att,
  treat_type = "post",
  estimators = ESTIMATORS
)

print(example_results %>% select(method, estimate, se, true_att, bias, p_value, se_fallback))
write_csv(example_results, "output/example_estimates.csv")

# =============================================================================
# PART 2: Power (vary delta)
# =============================================================================
power_results <- run_power_analysis(
  delta_values = c(0, 0.2, 0.4, 0.6, 0.8),
  n_sims = 200,
  estimators = ESTIMATORS,
  outcome_col = OUTCOME_COL,
  treat_type = "post",
  base_params = list(
    s = 0.9,
    effect_type = "constant",
    rho = 0.3,
    sigma_u = anchor$sigma_u,
    q = 1.0,
    decay_half_life = 3,
    mu_pos_cap = 30
  ),
  lat_range = LAT_RANGE,
  lon_range = LON_RANGE,
  start_date = START_DATE,
  end_date = END_DATE,
  treat_start = TREAT_START,
  near_cutoff_km = NEAR_CUTOFF,
  K = K_FACTORS
)

power_summary <- summarize_power_results(power_results)
write_csv(power_summary, "output/power_summary.csv")

# Report SE fallback rates
cat("\nSE Fallback Rates (power analysis):\n")
print(power_summary %>% select(method, delta, se_fallback_rate))

ggsave("figures/power_curves.pdf", plot_power_curves(power_summary), width = 7.2, height = 4.6)
ggsave("figures/bias_curve.pdf", plot_bias(power_summary), width = 7.2, height = 4.6)
ggsave("figures/coverage_curve.pdf", plot_coverage(power_summary), width = 7.2, height = 4.6)

cat("\nSaved: output/power_summary.csv and figures/*.pdf\n")

# =============================================================================
# PART 3: Scenarios
# =============================================================================
scenario_results <- run_scenario_analysis(
  n_sims = 100,
  estimators = ESTIMATORS,
  outcome_col = OUTCOME_COL,
  treat_type = "auto",
  lat_range = LAT_RANGE,
  lon_range = LON_RANGE,
  start_date = START_DATE,
  end_date = END_DATE,
  treat_start = TREAT_START,
  near_cutoff_km = NEAR_CUTOFF,
  K = K_FACTORS
)

scenario_summary <- summarize_scenario_results(scenario_results)
write_csv(scenario_summary, "output/scenario_summary.csv")

# Report SE fallback rates
cat("\nSE Fallback Rates (scenario analysis):\n")
print(scenario_summary %>% select(scenario, method, se_fallback_rate))

cat("\nSaved: output/scenario_summary.csv\n")
cat("DONE\n")
