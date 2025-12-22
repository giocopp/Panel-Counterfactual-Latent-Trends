#' ============================================================================
#' Monte Carlo Simulation (Power + Scenarios)
#' Clean comparison: TWFE vs Matrix Completion vs SynthDiD
#' Outcome: Y_any (binary indicator)
#' ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

USE_PARALLEL <- FALSE
if (requireNamespace("furrr", quietly = TRUE)) {
  suppressPackageStartupMessages(library(furrr))
  plan(multisession, workers = max(1, parallel::detectCores() - 1))
  USE_PARALLEL <- TRUE
  cat("Parallel processing enabled:", max(1, parallel::detectCores() - 1), "workers\n")
}

source("R/01_dgp.R")
source("R/02_calibration.R")
source("R/03_estimators.R")

AVAILABLE_ESTIMATORS <- c("twfe", "mc")
if (exists("SYNTHDID_AVAILABLE") && isTRUE(SYNTHDID_AVAILABLE)) {
  AVAILABLE_ESTIMATORS <- c(AVAILABLE_ESTIMATORS, "sdid")
}

ensure_anchor <- function(
  anchor_path = "data/dgp_anchor.rds",
  lat_range = c(31, 38),
  lon_range = c(10, 18),
  start_date = "2015-04-01",
  end_date = "2018-02-01",
  treat_start = "2017-05-01",
  near_cutoff_km = 200,
  K = 2,
  c_smooth = 0.5
) {
  if (file.exists(anchor_path)) {
    return(load_dgp_anchor(anchor_path))
  }

  cat("No anchor found at", anchor_path, "- building from IOM...\n")
  cal <- calibrate_from_iom(
    lat_range = lat_range, lon_range = lon_range,
    start_date = start_date, end_date = end_date,
    treat_start = treat_start, near_cutoff_km = near_cutoff_km,
    K = K, c_smooth = c_smooth
  )
  if (is.null(cal) || is.null(cal$anchor)) stop("Failed to build anchor.")
  cal$anchor
}

get_target_pre_rate <- function(
  grid,
  time_panel,
  iom_panel_path = "data/iom_grid_panel.rds",
  iom_path = NULL
) {
  if (file.exists(iom_panel_path)) {
    iom_panel <- readRDS(iom_panel_path)
    return(mean(iom_panel$Y_any[iom_panel$post_treat == 0], na.rm = TRUE))
  }
  if (is.null(iom_path)) iom_path <- find_iom_file()
  if (is.null(iom_path)) stop("Cannot compute target_pre_rate: IOM file not found.")

  iom <- read_iom_incidents(iom_path)
  iom_panel <- build_iom_grid_panel(iom, grid, time_panel)

  dir.create("data", showWarnings = FALSE, recursive = TRUE)
  saveRDS(iom_panel, iom_panel_path)

  mean(iom_panel$Y_any[iom_panel$post_treat == 0], na.rm = TRUE)
}

# =============================================================================
# SINGLE SIMULATION RUN
# =============================================================================
run_one_sim <- function(
  sim_id,
  scenario_params,
  estimators = AVAILABLE_ESTIMATORS,
  outcome_col = "Y_any",
  treat_type = c("auto", "post", "window"),
  trim_after_window = NULL,
  grid = NULL,
  time_panel = NULL,
  anchor = NULL
) {
  treat_type <- match.arg(treat_type)
  effect_type <- scenario_params$effect_type %||% "constant"

  treat_type_eff <- if (treat_type == "auto") {
    if (identical(effect_type, "window")) "window" else "post"
  } else {
    treat_type
  }
  if (is.null(trim_after_window)) trim_after_window <- (treat_type_eff == "window")

  data <- tryCatch(
    do.call(generate_panel_data, c(list(grid = grid, time_panel = time_panel, anchor = anchor, seed = sim_id), scenario_params)),
    error = function(e) {
      message("DGP failed for sim ", sim_id, ": ", e$message)
      NULL
    }
  )

  if (is.null(data)) {
    method_names <- c(twfe = "TWFE", mc = "Matrix Completion", sdid = "SynthDiD")
    return(tibble(
      sim_id = sim_id,
      method = method_names[estimators],
      estimate = NA_real_, se = NA_real_, p_value = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_, true_att = NA_real_,
      se_fallback = FALSE
    ))
  }

  true <- tryCatch(
    compute_true_estimands(data, outcome_col = outcome_col, treat_type = treat_type_eff, trim_after_window = trim_after_window),
    error = function(e) list(att = NA_real_)
  )
  true_att <- true$att

  out <- tibble()

  if ("twfe" %in% estimators) {
    res <- tryCatch(
      estimate_twfe(panel = data, outcome_col = outcome_col, treat_type = treat_type_eff, trim_after_window = trim_after_window),
      error = function(e) list(method = "TWFE", estimate = NA_real_, se = NA_real_, p_value = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, se_fallback = FALSE)
    )
    out <- bind_rows(out, tibble(
      sim_id = sim_id, method = res$method, estimate = res$estimate, se = res$se, p_value = res$p_value,
      ci_lower = res$ci_lower, ci_upper = res$ci_upper, true_att = true_att,
      se_fallback = res$se_fallback %||% FALSE
    ))
  }

  if ("mc" %in% estimators) {
    res <- tryCatch(
      estimate_mc(panel = data, outcome_col = outcome_col, treat_type = treat_type_eff, trim_after_window = trim_after_window),
      error = function(e) list(method = "Matrix Completion", estimate = NA_real_, se = NA_real_, p_value = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, se_fallback = FALSE)
    )
    out <- bind_rows(out, tibble(
      sim_id = sim_id, method = res$method, estimate = res$estimate, se = res$se, p_value = res$p_value,
      ci_lower = res$ci_lower, ci_upper = res$ci_upper, true_att = true_att,
      se_fallback = res$se_fallback %||% FALSE
    ))
  }

  if ("sdid" %in% estimators) {
    res <- tryCatch(
      estimate_sdid(panel = data, outcome_col = outcome_col, treat_type = treat_type_eff, trim_after_window = trim_after_window),
      error = function(e) list(method = "SynthDiD", estimate = NA_real_, se = NA_real_, p_value = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, se_fallback = FALSE)
    )
    out <- bind_rows(out, tibble(
      sim_id = sim_id, method = res$method, estimate = res$estimate, se = res$se, p_value = res$p_value,
      ci_lower = res$ci_lower, ci_upper = res$ci_upper, true_att = true_att,
      se_fallback = res$se_fallback %||% FALSE
    ))
  }

  out
}

# =============================================================================
# POWER ANALYSIS (vary delta)
# =============================================================================
run_power_analysis <- function(
  delta_values = c(0, 0.2, 0.4, 0.6, 0.8),
  n_sims = 200,
  estimators = AVAILABLE_ESTIMATORS,
  outcome_col = "Y_any",
  treat_type = c("auto", "post", "window"),
  trim_after_window = NULL,
  base_params = list(
    s = 0.9,
    effect_type = "constant",
    rho = 0.3,
    sigma_u = NULL,
    q = 1.0,
    decay_half_life = 3,
    mu_pos_cap = 30
  ),
  lat_range = c(31, 38),
  lon_range = c(10, 18),
  start_date = "2015-04-01",
  end_date = "2018-02-01",
  treat_start = "2017-05-01",
  near_cutoff_km = 200,
  K = 2
) {
  treat_type <- match.arg(treat_type)

  grid <- generate_grid(lat_range = lat_range, lon_range = lon_range, near_cutoff_km = near_cutoff_km)
  time_panel <- generate_time_panel(start_date = start_date, end_date = end_date, treat_start = treat_start)
  anchor <- ensure_anchor(
    lat_range = lat_range, lon_range = lon_range, start_date = start_date, end_date = end_date,
    treat_start = treat_start, near_cutoff_km = near_cutoff_km, K = K
  )

  target_pre_rate <- get_target_pre_rate(grid = grid, time_panel = time_panel)
  cat("Target pre-rate (IOM, Y_any):", signif(target_pre_rate, 4), "\n")

  if (is.null(base_params$auto_level_shift)) base_params$auto_level_shift <- TRUE
  if (is.null(base_params$target_pre_rate)) base_params$target_pre_rate <- target_pre_rate

  sim_grid <- expand.grid(sim_id = 1:n_sims, delta = delta_values) %>% as_tibble()
  cat("Running", nrow(sim_grid), "simulations\n")
  cat("Outcome:", outcome_col, "\n")
  cat("Estimators:", paste(estimators, collapse = ", "), "\n")
  cat("Parallel:", USE_PARALLEL, "\n\n")

  runner <- function(sim_id, delta) {
    params <- c(list(delta = delta), base_params)
    run_one_sim(
      sim_id = sim_id, scenario_params = params, estimators = estimators,
      outcome_col = outcome_col, treat_type = treat_type, trim_after_window = trim_after_window,
      grid = grid, time_panel = time_panel, anchor = anchor
    ) %>%
      mutate(delta = delta)
  }

  if (USE_PARALLEL) {
    furrr::future_pmap_dfr(
      sim_grid,
      runner,
      .options = furrr::furrr_options(seed = TRUE, packages = c("fixest", "tidyverse", "synthdid", "fect")),
      .progress = TRUE
    )
  } else {
    purrr::pmap_dfr(sim_grid, runner)
  }
}

# =============================================================================
# SCENARIO ANALYSIS
# =============================================================================
run_scenario_analysis <- function(
  scenarios = NULL,
  n_sims = 200,
  estimators = AVAILABLE_ESTIMATORS,
  outcome_col = "Y_any",
  treat_type = c("auto", "post", "window"),
  trim_after_window = NULL,
  lat_range = c(31, 38),
  lon_range = c(10, 18),
  start_date = "2015-04-01",
  end_date = "2018-02-01",
  treat_start = "2017-05-01",
  near_cutoff_km = 200,
  K = 2
) {
  if (is.null(scenarios)) scenarios <- dgp_scenarios()
  treat_type <- match.arg(treat_type)

  grid <- generate_grid(lat_range = lat_range, lon_range = lon_range, near_cutoff_km = near_cutoff_km)
  time_panel <- generate_time_panel(start_date = start_date, end_date = end_date, treat_start = treat_start)
  anchor <- ensure_anchor(
    lat_range = lat_range, lon_range = lon_range, start_date = start_date, end_date = end_date,
    treat_start = treat_start, near_cutoff_km = near_cutoff_km, K = K
  )

  map_dfr(names(scenarios), function(scenario_name) {
    cat("\n=== Scenario:", scenario_name, "===\n")
    params <- scenarios[[scenario_name]]

    runner <- function(sim_id) {
      run_one_sim(
        sim_id = sim_id, scenario_params = params, estimators = estimators,
        outcome_col = outcome_col, treat_type = treat_type, trim_after_window = trim_after_window,
        grid = grid, time_panel = time_panel, anchor = anchor
      )
    }

    sim_results <- if (USE_PARALLEL) {
      furrr::future_map_dfr(
        1:n_sims,
        runner,
        .options = furrr::furrr_options(seed = TRUE, packages = c("fixest", "tidyverse", "synthdid", "fect")),
        .progress = TRUE
      )
    } else {
      purrr::map_dfr(1:n_sims, runner)
    }

    sim_results %>% mutate(scenario = scenario_name)
  })
}

# =============================================================================
# SUMMARY + PLOTTING
# =============================================================================
summarize_power_results <- function(results, alpha = 0.05) {
  results %>%
    filter(!is.na(estimate)) %>%
    group_by(method, delta) %>%
    summarise(
      n_sims = n(),
      true_att = mean(true_att, na.rm = TRUE),
      mean_estimate = mean(estimate, na.rm = TRUE),
      bias = mean(estimate - true_att, na.rm = TRUE),
      rmse = sqrt(mean((estimate - true_att)^2, na.rm = TRUE)),
      coverage = mean(ci_lower <= true_att & ci_upper >= true_att, na.rm = TRUE),
      rejection_rate = mean(p_value < alpha, na.rm = TRUE),
      se_fallback_rate = mean(se_fallback, na.rm = TRUE),
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
      se_fallback_rate = mean(se_fallback, na.rm = TRUE),
      .groups = "drop"
    )
}

METHOD_COLORS <- c(
  "TWFE" = "#E69F00",
  "Matrix Completion" = "#56B4E9",
  "SynthDiD" = "#CC79A7"
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
    labs(title = "Power curves", x = "True ATT", y = "Rejection rate", color = "Method", shape = "Method") +
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
    labs(title = "Bias by true ATT", x = "True ATT", y = "Bias", color = "Method", shape = "Method") +
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
    labs(title = "95% CI coverage", x = "True ATT", y = "Coverage", color = "Method", shape = "Method") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}
