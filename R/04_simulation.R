#' ============================================================================
#' Monte Carlo Simulation for Power Analysis
#'
#' Focus: panel counterfactual estimators under latent trends (interactive FE)
#' Outcome: default Y_any (deadly-event indicator)
#' Treatment: D_it = near_libya_i × post_treat_t  (post_treat = May 2017)
#'
#' Metrics:
#'   - rejection rate (power/size)
#'   - bias, RMSE
#'   - coverage (nominal 95% CI)
#' ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(lubridate)
  library(scales)
})

source("R/01_dgp.R")
source("R/03_estimators.R")

run_one_sim <- function(
    sim_id,
    scenario_params,
    estimators = c("twfe", "mc", "sdid", "trop"),
    outcome_col = "Y_any",
    grid = NULL,
    time_panel = NULL
) {
  data <- do.call(generate_panel_data, c(
    list(grid = grid, time_panel = time_panel, seed = sim_id),
    scenario_params
  ))

  true <- compute_true_estimands(data, outcome_col = outcome_col)
  true_att <- true$att

  results <- tibble()

  if ("twfe" %in% estimators) {
    res <- tryCatch(
      estimate_twfe(data, outcome_col = outcome_col),
      error = function(e) list(estimate = NA_real_, se = NA_real_, p_value = NA_real_,
                               ci_lower = NA_real_, ci_upper = NA_real_)
    )
    results <- bind_rows(results, tibble(
      sim_id = sim_id,
      method = "TWFE",
      estimate = res$estimate,
      se = res$se,
      p_value = res$p_value,
      ci_lower = res$ci_lower,
      ci_upper = res$ci_upper,
      true_att = true_att
    ))
  }

  if ("mc" %in% estimators) {
    res <- tryCatch(
      estimate_mc(data, outcome_col = outcome_col),
      error = function(e) list(estimate = NA_real_, se = NA_real_, p_value = NA_real_,
                               ci_lower = NA_real_, ci_upper = NA_real_)
    )
    results <- bind_rows(results, tibble(
      sim_id = sim_id,
      method = "Matrix Completion",
      estimate = res$estimate,
      se = res$se,
      p_value = res$p_value,
      ci_lower = res$ci_lower,
      ci_upper = res$ci_upper,
      true_att = true_att
    ))
  }

  if ("sdid" %in% estimators) {
    res <- tryCatch(
      estimate_sdid(data, outcome_col = outcome_col),
      error = function(e) list(estimate = NA_real_, se = NA_real_, p_value = NA_real_,
                               ci_lower = NA_real_, ci_upper = NA_real_)
    )
    results <- bind_rows(results, tibble(
      sim_id = sim_id,
      method = "Synthetic DiD",
      estimate = res$estimate,
      se = res$se,
      p_value = res$p_value,
      ci_lower = res$ci_lower,
      ci_upper = res$ci_upper,
      true_att = true_att
    ))
  }

  if ("trop" %in% estimators) {
    res <- tryCatch(
      estimate_trop(data, outcome_col = outcome_col),
      error = function(e) list(estimate = NA_real_, se = NA_real_, p_value = NA_real_,
                               ci_lower = NA_real_, ci_upper = NA_real_)
    )
    results <- bind_rows(results, tibble(
      sim_id = sim_id,
      method = "TROP",
      estimate = res$estimate,
      se = res$se,
      p_value = res$p_value,
      ci_lower = res$ci_lower,
      ci_upper = res$ci_upper,
      true_att = true_att
    ))
  }

  results
}

run_power_analysis <- function(
    delta_values = c(0, 0.2, 0.4, 0.6, 0.8),
    n_sims = 200,
    estimators = c("twfe", "mc", "sdid", "trop"),
    outcome_col = "Y_any",
    base_params = list(
      factor_strength = 0.5,
      underreporting_rate = 0.0,
      short_lived = TRUE
    ),
    baseline_mortality = 0.02,
    seasonality_index = NULL
) {

  grid <- generate_grid(near_cutoff_km = 200)
  time_panel <- generate_time_panel(seasonality_index = seasonality_index)

  sim_grid <- expand.grid(
    sim_id = 1:n_sims,
    delta = delta_values
  ) %>%
    as_tibble()

  cat("Running", nrow(sim_grid), "simulations\n")
  cat("Outcome:", outcome_col, "\n\n")

  pmap_dfr(sim_grid, function(sim_id, delta) {
    params <- c(list(delta = delta, baseline_mortality = baseline_mortality), base_params)
    run_one_sim(sim_id, params, estimators = estimators, outcome_col = outcome_col,
                grid = grid, time_panel = time_panel) %>%
      mutate(delta = delta)
  })
}

run_scenario_analysis <- function(
    scenarios = NULL,
    n_sims = 200,
    estimators = c("twfe", "mc", "sdid", "trop"),
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

    sim_results <- map_dfr(1:n_sims, function(sim_id) {
      if (sim_id %% 50 == 0) cat("  Sim", sim_id, "\n")
      run_one_sim(sim_id, params, estimators = estimators, outcome_col = outcome_col,
                  grid = grid, time_panel = time_panel)
    })

    sim_results %>% mutate(scenario = scenario_name)
  })
}

summarize_power_results <- function(results, alpha = 0.05) {
  results %>%
    group_by(method, delta) %>%
    summarise(
      n_sims = dplyr::n(),
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
    group_by(scenario, method) %>%
    summarise(
      n_sims = dplyr::n(),
      true_att = mean(true_att, na.rm = TRUE),
      bias = mean(estimate - true_att, na.rm = TRUE),
      rmse = sqrt(mean((estimate - true_att)^2, na.rm = TRUE)),
      coverage = mean(ci_lower <= true_att & ci_upper >= true_att, na.rm = TRUE),
      power = mean(p_value < alpha, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_power_curves <- function(summary) {
  summary %>%
    ggplot(aes(x = true_att, y = rejection_rate, color = method, shape = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.8) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0.80, linetype = "dotted", color = "gray40") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title = "Power curves (x-axis is implied true ATT)",
      subtitle = "Dashed: 5% size. Dotted: 80% power target.",
      x = "True ATT (E[Y(1)-Y(0)] among treated post)",
      y = "Rejection rate",
      color = "Method",
      shape = "Method"
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
    labs(
      title = "Bias by implied true ATT",
      x = "True ATT",
      y = "Bias (estimate - true)",
      color = "Method",
      shape = "Method"
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
    scale_y_continuous(limits = c(0.7, 1), labels = scales::percent) +
    labs(
      title = "95% CI coverage by implied true ATT",
      x = "True ATT",
      y = "Coverage",
      color = "Method",
      shape = "Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}
