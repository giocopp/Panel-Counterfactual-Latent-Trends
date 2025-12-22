#' ============================================================================
#' Calibration using IOM MM data (for Central Mediterranean)
#'
#' Purpose: set DGP baselines to resemble the sparsity + seasonality observed in
#'          IOM Missing Migrants data aggregated to the same grid×month panel.
#'
#' This is intentionally lightweight calibration for a methods project:
#'   - Target 1: pre-period mean of Y_any = 1{deaths > 0} across cell-months
#'   - Target 2: month-of-year multipliers for event probability (seasonality)
#' ============================================================================

library(tidyverse)
library(lubridate)

`%||%` <- function(x, y) if (is.null(x)) y else x

find_iom_file <- function() {
  candidates <- c(
    "data/iom_mediterranean_2014_2020.csv",
    "iom_mediterranean_2014_2020.csv"
  )
  for (p in candidates) {
    if (file.exists(p)) {
      message("Using IOM data: ", p)
      return(p)
    }
  }
  return(NULL)
}

read_iom_incidents <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(
      incident_date = as.Date(`Incident date`),
      deaths = readr::parse_number(as.character(`No. dead/missing`)),
      lat = as.numeric(Latitude),
      lon = as.numeric(Longitude)
    ) %>%
    filter(
      !is.na(incident_date),
      !is.na(lat), !is.na(lon),
      !is.na(deaths),
      str_detect(Route, "Central Mediterranean")
    ) %>%
    mutate(month_date = floor_date(incident_date, unit = "month"))
}

# Build a full grid×month panel on the SAME grid definition as the DGP
build_iom_grid_panel <- function(iom_incidents, grid, time_panel) {
  # Map incidents to 1° cells (by floor of coordinates)
  iom_cells <- iom_incidents %>%
    mutate(
      lat_cell = floor(lat),
      lon_cell = floor(lon),
      unit_key = paste0(lat_cell, "_", lon_cell)
    )

  # Create a key for DGP cells (same convention)
  grid_key <- grid %>%
    mutate(unit_key = paste0(floor(lat_center), "_", floor(lon_center))) %>%
    select(unit_id, unit_key)

  # Keep incidents that land inside the DGP grid keys
  iom_cells <- iom_cells %>%
    inner_join(grid_key, by = "unit_key")

  # Aggregate to cell×month
  cell_month <- iom_cells %>%
    group_by(unit_id, month_date) %>%
    summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop")

  # Full balanced panel
  panel <- tidyr::crossing(
    grid %>% select(unit_id),
    time_panel %>% select(date, time_id, month, post_treat)
  ) %>%
    rename(month_date = date) %>%
    left_join(cell_month, by = c("unit_id", "month_date")) %>%
    mutate(
      deaths = replace_na(deaths, 0),
      Y_any = as.integer(deaths > 0),
      Y_ihs = asinh(deaths)
    )

  panel
}

compute_calibration_targets <- function(iom_panel, treat_start = "2017-05-01") {
  treat_start <- as.Date(treat_start)

  pre <- iom_panel %>%
    filter(month_date < treat_start)

  target_event_rate <- mean(pre$Y_any)

  seasonality_index <- pre %>%
    group_by(month) %>%
    summarise(rate = mean(Y_any), .groups = "drop") %>%
    mutate(index = rate / mean(rate)) %>%
    arrange(month) %>%
    pull(index)

  list(
    target_event_rate = target_event_rate,
    seasonality_index = seasonality_index
  )
}

# Calibrate baseline_mortality so that simulated pre-period mean(Y_any) matches IOM target.
# We do this with a small uniroot search on baseline_mortality.
calibrate_baseline_mortality <- function(
  target_event_rate,
  grid,
  time_panel,
  baseline_attempts = 60,
  factor_strength = 0.5,
  near_cutoff_km = 200,
  seed = 123
) {
  obj <- function(bm) {
    sim <- generate_panel_data(
      grid = grid,
      time_panel = time_panel,
      delta = 0.0, # no treatment effect during calibration
      baseline_attempts = baseline_attempts,
      baseline_mortality = bm,
      factor_strength = factor_strength,
      underreporting_rate = 0.0,
      coord_missing_rate = 0.0,
      short_lived = FALSE,
      near_cutoff_km = near_cutoff_km,
      seed = seed
    )

    pre <- sim %>% filter(post_treat == 0)
    mean(pre$Y_any) - target_event_rate
  }

  # Search bounds: keep small; event rates explode otherwise
  lower <- 1e-7
  upper <- 0.01

  f_lower <- obj(lower)
  f_upper <- obj(upper)

  if (is.na(f_lower) || is.na(f_upper) || f_lower * f_upper > 0) {
    # Fall back to a simple grid search if uniroot can't bracket
    candidates <- 10^seq(-7, -2, length.out = 25)
    vals <- purrr::map_dbl(candidates, obj)
    bm_best <- candidates[which.min(abs(vals))]
    return(list(baseline_mortality = bm_best, method = "grid_search"))
  }

  bm <- uniroot(obj, interval = c(lower, upper), tol = 1e-6)$root
  list(baseline_mortality = bm, method = "uniroot")
}

# Main entry point used by run_analysis.R
calibrate_from_iom <- function(
  iom_path = NULL,
  lat_range = c(31, 38),
  lon_range = c(10, 18),
  start_date = "2015-04-01",
  end_date = "2018-02-01",
  treat_start = "2017-05-01",
  baseline_attempts = 60,
  factor_strength = 0.5,
  near_cutoff_km = 200
) {
  if (is.null(iom_path)) iom_path <- find_iom_file()
  if (is.null(iom_path)) {
    message("No IOM incident file found; skipping calibration.")
    return(NULL)
  }

  iom <- read_iom_incidents(iom_path)

  # Build DGP grid + time panel for mapping and seasonality extraction
  grid <- generate_grid(lat_range = lat_range, lon_range = lon_range, near_cutoff_km = near_cutoff_km)

  # Temporary time panel without seasonality_index (we'll compute it from IOM)
  time_panel_raw <- generate_time_panel(
    start_date = start_date,
    end_date = end_date,
    treat_start = treat_start
  )

  iom_panel <- build_iom_grid_panel(iom, grid, time_panel_raw)
  targets <- compute_calibration_targets(iom_panel, treat_start = treat_start)

  # Use IOM seasonality_index to build the time panel used in simulation
  time_panel <- generate_time_panel(
    start_date = start_date,
    end_date = end_date,
    treat_start = treat_start,
    seasonality_index = targets$seasonality_index
  )

  bm <- calibrate_baseline_mortality(
    target_event_rate = targets$target_event_rate,
    grid = grid,
    time_panel = time_panel,
    baseline_attempts = baseline_attempts,
    factor_strength = factor_strength,
    near_cutoff_km = near_cutoff_km
  )

  result <- list(
    iom_path = iom_path,
    target_event_rate = targets$target_event_rate,
    seasonality_index = targets$seasonality_index,
    baseline_mortality = bm$baseline_mortality,
    calibration_method = bm$method
  )

  saveRDS(result, "data/calibration_targets.rds")
  message("Calibration saved to data/calibration_targets.rds")

  result
}

# =============================================================================
# Diagnostics helpers (figures for the writeup)
# =============================================================================

#' Create simple calibration diagnostic figures.
#'
#' Produces:
#'   - figures/calibration_seasonality.pdf
#'   - figures/calibration_pre_rate_check.pdf
#'   - figures/calibration_iom_monthly_rate.pdf (if IOM file is available)
make_calibration_figures <- function(
  calibration,
  grid,
  time_panel,
  treat_start = "2017-05-01",
  fig_dir = "figures",
  n_rep = 30,
  seed = 202
) {
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  treat_start <- as.Date(treat_start)

  # ---- 1) Seasonality multipliers (month-of-year)
  if (!is.null(calibration$seasonality_index)) {
    df_seas <- tibble(
      month = 1:12,
      index = as.numeric(calibration$seasonality_index),
      month_lab = factor(month.abb[month], levels = month.abb)
    )
    p_seas <- ggplot(df_seas, aes(x = month_lab, y = index)) +
      geom_col() +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(
        title = "Calibration: seasonality multipliers (pre-period)",
        x = NULL,
        y = "Multiplier (relative to mean)"
      ) +
      theme_minimal()
    ggsave(
      filename = file.path(fig_dir, "calibration_seasonality.pdf"),
      plot = p_seas,
      width = 7,
      height = 3.8
    )
  }

  # ---- 2) Pre-period event-rate sanity check under the calibrated baseline
  # Run small repeated simulations with delta = 0 and compare mean(Y_any) pre.
  set.seed(seed)
  pre_means <- rep(NA_real_, n_rep)
  for (r in seq_len(n_rep)) {
    sim <- generate_panel_data(
      grid = grid,
      time_panel = time_panel,
      delta = 0,
      baseline_mortality = calibration$baseline_mortality,
      short_lived = FALSE,
      seed = seed + r
    )
    pre_means[r] <- mean(sim$Y_any[sim$date < treat_start], na.rm = TRUE)
  }

  df_check <- tibble(pre_mean = pre_means)
  target <- calibration$target_event_rate %||% mean(pre_means, na.rm = TRUE)
  p_check <- ggplot(df_check, aes(x = pre_mean)) +
    geom_histogram(bins = 20) +
    geom_vline(xintercept = target, linetype = "dashed") +
    labs(
      title = "Calibration: simulated pre-period event rate (delta = 0)",
      subtitle = paste0("Dashed line = target ", signif(target, 3)),
      x = "Mean(Y_any) in pre-period",
      y = "Count"
    ) +
    theme_minimal()
  ggsave(
    filename = file.path(fig_dir, "calibration_pre_rate_check.pdf"),
    plot = p_check,
    width = 7,
    height = 3.8
  )

  # ---- 3) If IOM raw file is present, plot monthly pre-period event rates
  iom_path <- calibration$iom_path %||% find_iom_file()
  if (!is.null(iom_path) && file.exists(iom_path)) {
    iom <- read_iom_incidents(iom_path)
    iom_panel <- build_iom_grid_panel(iom, grid, time_panel)
    df_iom <- iom_panel %>%
      filter(month_date < treat_start) %>%
      group_by(month) %>%
      summarise(rate = mean(Y_any), .groups = "drop") %>%
      mutate(month_lab = factor(month.abb[month], levels = month.abb))

    p_iom <- ggplot(df_iom, aes(x = month_lab, y = rate)) +
      geom_col() +
      labs(
        title = "IOM (pre-period): monthly deadly-event probability",
        x = NULL,
        y = "Mean(Y_any)"
      ) +
      theme_minimal()

    ggsave(
      filename = file.path(fig_dir, "calibration_iom_monthly_rate.pdf"),
      plot = p_iom,
      width = 7,
      height = 3.8
    )
  }

  invisible(TRUE)
}
