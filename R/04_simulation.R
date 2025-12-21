#' ============================================================================
#' Monte Carlo Simulation for Power Analysis
#'
#' Compares TWFE, Matrix Completion, TROP, and SynthDiD
#' estimators under various DGP configurations with latent factor structure.
#' ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
})

# Enable parallel processing if available
USE_PARALLEL <- FALSE
if (requireNamespace("furrr", quietly = TRUE)) {
  library(furrr)
  plan(multisession, workers = max(1, parallel::detectCores() - 1))
  USE_PARALLEL <- TRUE
  cat("Parallel processing enabled:", parallel::detectCores() - 1, "workers\n")
}

source("R/01_dgp.R")
source("R/03_estimators.R")

# Check which estimators are available
AVAILABLE_ESTIMATORS <- c("twfe", "mc", "trop")
if (exists("SYNTHDID_AVAILABLE") && SYNTHDID_AVAILABLE) {
  AVAILABLE_ESTIMATORS <- c(AVAILABLE_ESTIMATORS, "sdid")
}

# =============================================================================
# SINGLE SIMULATION RUN
# =============================================================================

run_one_sim <- function(
  sim_id,
  scenario_params,
  estimators = AVAILABLE_ESTIMATORS,
  outcome_col = "Y_any",
  grid = NULL,
  time_panel = NULL
) {
  # Generate data using DGP
  data <- tryCatch(
    {
      do.call(generate_panel_data, c(
        list(grid = grid, time_panel = time_panel, seed = sim_id),
        scenario_params
      ))
    },
    error = function(e) {
      message("DGP failed for sim ", sim_id, ": ", e$message)
      return(NULL)
    }
  )

  if (is.null(data)) {
    method_names <- c(
      twfe = "TWFE",
      mc = "Matrix Completion",
      trop = "TROP",
      sdid = "SynthDiD"
    )
    return(tibble(
      sim_id = sim_id,
      method = method_names[estimators],
      estimate = NA_real_,
      se = NA_real_,
      p_value = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      true_att = NA_real_
    ))
  }

  # Compute true ATT
  true <- tryCatch(
    compute_true_estimands(data, outcome_col = outcome_col),
    error = function(e) list(att = NA_real_)
  )
  true_att <- true$att

  results <- tibble()

  if ("twfe" %in% estimators) {
    res <- tryCatch(
      estimate_twfe(data, outcome_col = outcome_col),
      error = function(e) {
        list(
          estimate = NA_real_, se = NA_real_, p_value = NA_real_,
          ci_lower = NA_real_, ci_upper = NA_real_
        )
      }
    )
    results <- bind_rows(results, tibble(
      sim_id = sim_id, method = "TWFE", estimate = res$estimate, se = res$se,
      p_value = res$p_value, ci_lower = res$ci_lower, ci_upper = res$ci_upper,
      true_att = true_att
    ))
  }

  if ("mc" %in% estimators) {
    res <- tryCatch(
      estimate_mc(data, outcome_col = outcome_col),
      error = function(e) {
        list(
          estimate = NA_real_, se = NA_real_, p_value = NA_real_,
          ci_lower = NA_real_, ci_upper = NA_real_
        )
      }
    )
    results <- bind_rows(results, tibble(
      sim_id = sim_id, method = "Matrix Completion", estimate = res$estimate, se = res$se,
      p_value = res$p_value, ci_lower = res$ci_lower, ci_upper = res$ci_upper,
      true_att = true_att
    ))
  }

  if ("trop" %in% estimators) {
    res <- tryCatch(
      estimate_trop(data, outcome_col = outcome_col),
      error = function(e) {
        list(
          estimate = NA_real_, se = NA_real_, p_value = NA_real_,
          ci_lower = NA_real_, ci_upper = NA_real_
        )
      }
    )
    results <- bind_rows(results, tibble(
      sim_id = sim_id, method = "TROP", estimate = res$estimate, se = res$se,
      p_value = res$p_value, ci_lower = res$ci_lower, ci_upper = res$ci_upper,
      true_att = true_att
    ))
  }

  if ("sdid" %in% estimators) {
    res <- tryCatch(
      estimate_sdid(data, outcome_col = outcome_col),
      error = function(e) {
        list(
          estimate = NA_real_, se = NA_real_, p_value = NA_real_,
          ci_lower = NA_real_, ci_upper = NA_real_
        )
      }
    )
    results <- bind_rows(results, tibble(
      sim_id = sim_id, method = "SynthDiD", estimate = res$estimate, se = res$se,
      p_value = res$p_value, ci_lower = res$ci_lower, ci_upper = res$ci_upper,
      true_att = true_att
    ))
  }

  results
}

# =============================================================================
# POWER ANALYSIS
# =============================================================================

run_power_analysis <- function(
  delta_values = c(0, 0.2, 0.4, 0.6, 0.8),
  n_sims = 200,
  estimators = AVAILABLE_ESTIMATORS,
  outcome_col = "Y_any",
  base_params = list(
    factor_strength = 0.5,
    loading_correlation = 0.5,  # Correlated loadings -> violates parallel trends
    underreporting_rate = 0.0,
    short_lived = TRUE
  ),
  baseline_mortality = 0.02,
  seasonality_index = NULL
) {
  # Generate spatial grid and time panel once
  grid <- generate_grid(near_cutoff_km = 200)
  time_panel <- generate_time_panel(seasonality_index = seasonality_index)

  # Create simulation grid
  sim_grid <- expand.grid(
    sim_id = 1:n_sims,
    delta = delta_values
  ) %>% as_tibble()

  cat("Running", nrow(sim_grid), "simulations\n")
  cat("Outcome:", outcome_col, "\n")
  cat("Estimators:", paste(estimators, collapse = ", "), "\n")
  cat("Parallel:", USE_PARALLEL, "\n\n")

  # Choose mapping function based on parallel availability
  if (USE_PARALLEL) {
    results <- furrr::future_pmap_dfr(
      sim_grid,
      function(sim_id, delta) {
        params <- c(list(delta = delta, baseline_mortality = baseline_mortality), base_params)
        run_one_sim(sim_id, params,
          estimators = estimators, outcome_col = outcome_col,
          grid = grid, time_panel = time_panel
        ) %>%
          mutate(delta = delta)
      },
      .options = furrr_options(seed = TRUE, packages = c("fixest", "tidyverse", "synthdid")),
      .progress = TRUE
    )
  } else {
    results <- purrr::pmap_dfr(
      sim_grid,
      function(sim_id, delta) {
        params <- c(list(delta = delta, baseline_mortality = baseline_mortality), base_params)
        run_one_sim(sim_id, params,
          estimators = estimators, outcome_col = outcome_col,
          grid = grid, time_panel = time_panel
        ) %>%
          mutate(delta = delta)
      },
      .progress = TRUE
    )
  }

  results
}

# =============================================================================
# SCENARIO ANALYSIS
# =============================================================================

run_scenario_analysis <- function(
  scenarios = NULL,
  n_sims = 200,
  estimators = AVAILABLE_ESTIMATORS,
  outcome_col = "Y_any",
  baseline_mortality = 0.02,
  seasonality_index = NULL
) {
  if (is.null(scenarios)) scenarios <- dgp_scenarios()

  grid <- generate_grid(near_cutoff_km = 200)
  time_panel <- generate_time_panel(seasonality_index = seasonality_index)

  map_dfr(names(scenarios), function(scenario_name) {
    cat("\n=== Scenario:", scenario_name, "===\n")
    params <- scenarios[[scenario_name]]
    params$baseline_mortality <- baseline_mortality

    if (USE_PARALLEL) {
      sim_results <- furrr::future_map_dfr(
        1:n_sims,
        function(sim_id) {
          run_one_sim(sim_id, params,
            estimators = estimators, outcome_col = outcome_col,
            grid = grid, time_panel = time_panel
          )
        },
        .options = furrr_options(seed = TRUE, packages = c("fixest", "tidyverse", "synthdid")),
        .progress = TRUE
      )
    } else {
      sim_results <- purrr::map_dfr(
        1:n_sims,
        function(sim_id) {
          run_one_sim(sim_id, params,
            estimators = estimators, outcome_col = outcome_col,
            grid = grid, time_panel = time_panel
          )
        },
        .progress = TRUE
      )
    }

    sim_results %>% mutate(scenario = scenario_name)
  })
}

# =============================================================================
# SUMMARY FUNCTIONS
# =============================================================================

summarize_power_results <- function(results, alpha = 0.05) {
  results %>%
    filter(!is.na(estimate)) %>%
    group_by(method, delta) %>%
    summarise(
      n_sims = n(),
      n_valid = sum(!is.na(estimate)),
      true_att = mean(true_att, na.rm = TRUE),
      mean_estimate = mean(estimate, na.rm = TRUE),
      bias = mean(estimate - true_att, na.rm = TRUE),
      rmse = sqrt(mean((estimate - true_att)^2, na.rm = TRUE)),
      coverage = mean(ci_lower <= true_att & ci_upper >= true_att, na.rm = TRUE),
      rejection_rate = mean(p_value < alpha, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(test_type = ifelse(abs(true_att) < 1e-8, "Size", "Power"))
}

summarize_scenario_results <- function(results, alpha = 0.05) {
  results %>%
    filter(!is.na(estimate)) %>%
    group_by(scenario, method) %>%
    summarise(
      n_sims = n(),
      true_att = mean(true_att, na.rm = TRUE),
      bias = mean(estimate - true_att, na.rm = TRUE),
      rmse = sqrt(mean((estimate - true_att)^2, na.rm = TRUE)),
      coverage = mean(ci_lower <= true_att & ci_upper >= true_att, na.rm = TRUE),
      power = mean(p_value < alpha, na.rm = TRUE),
      .groups = "drop"
    )
}

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

# Color palette for 4 methods (colorblind-friendly)
METHOD_COLORS <- c(
  "TWFE" = "#E69F00",           # Orange
  "Matrix Completion" = "#56B4E9", # Sky blue
  "TROP" = "#009E73",           # Teal
  "SynthDiD" = "#CC79A7"        # Pink
)

plot_power_curves <- function(summary) {
  summary %>%
    ggplot(aes(x = true_att, y = rejection_rate, color = method, shape = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.8) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0.80, linetype = "dotted", color = "gray40") +
    scale_y_continuous(limits = c(0, 1), labels = percent) +
    scale_color_manual(values = METHOD_COLORS) +
    labs(
      title = "Power curves",
      x = "True ATT",
      y = "Rejection rate",
      color = "Method", shape = "Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

plot_bias <- function(summary) {
  summary %>%
    ggplot(aes(x = true_att, y = bias, color = method, shape = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = METHOD_COLORS) +
    labs(
      title = "Bias by true ATT",
      x = "True ATT",
      y = "Bias",
      color = "Method", shape = "Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

plot_coverage <- function(summary) {
  summary %>%
    ggplot(aes(x = true_att, y = coverage, color = method, shape = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.8) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray40") +
    scale_y_continuous(limits = c(0.5, 1), labels = percent) +
    scale_color_manual(values = METHOD_COLORS) +
    labs(
      title = "95% CI coverage",
      x = "True ATT",
      y = "Coverage",
      color = "Method", shape = "Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

plot_scenario_comparison <- function(summary) {
  summary %>%
    select(scenario, method, rmse, coverage, power) %>%
    pivot_longer(cols = c(rmse, coverage, power), names_to = "metric") %>%
    ggplot(aes(x = method, y = value, fill = scenario)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~metric, scales = "free_y") +
    labs(
      title = "Estimator performance by scenario",
      x = "Method",
      y = "Value",
      fill = "Scenario"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}
