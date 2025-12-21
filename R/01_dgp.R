#' ============================================================================
#' Data Generating Process (DGP): Spatial Panel inspired by Zambiasi-Albarosa
#'
#' Creates a balanced unit×month panel inspired by spatial analyses of fatal
#' incidents on the Central Mediterranean route.
#'
#' Units: 1°×1° grid cells (cell centroids)
#' Time: monthly
#'
#' Key design features for causal-ML power analysis:
#'   - two-way fixed effects + interactive fixed effects (latent trends)
#'   - sparse outcomes (rare deadly events)
#'   - underreporting / measurement error
#'   - policy shock with stronger exposure for Libya-proximate cells
#'
#' IMPORTANT: The `loading_correlation` parameter controls whether factor
#' loadings are correlated with treatment assignment (near_libya). When > 0,
#' treated units have systematically different loadings, violating parallel
#' trends and biasing TWFE. This is the key mechanism that makes robust
#' estimators (MC, TROP) valuable.
#'
#' Primary outcome for estimator comparisons:
#'   Y_any = 1{ deaths_observed > 0 }  (cell-month deadly-event indicator)
#' Secondary outcome (robustness):
#'   Y_ihs = asinh(deaths_observed)
#' ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
})

# =============================================================================
# TIME PANEL
# =============================================================================

generate_time_panel <- function(
  start_date = "2015-04-01",
  end_date = "2018-02-01",
  mou_date = "2017-02-01",
  treat_start = "2017-05-01", # operational onset (patrol boats)
  other_policy = "2017-07-01" # e.g., NGO code of conduct
  , seasonality_index = NULL # optional 12-length vector (month multipliers)
) {
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "month")
  mou_date <- as.Date(mou_date)
  treat_start <- as.Date(treat_start)
  other_policy <- as.Date(other_policy)

  tibble(
    date = dates,
    year = year(dates),
    month = month(dates),
    time_id = seq_along(dates), # avoid dplyr::n()/row_number() outside verbs
    # Policy indicators
    post_mou = as.integer(dates >= mou_date),
    post_treat = as.integer(dates >= treat_start),
    # Short-lived "effect window" (May–June 2017)
    effect_window = as.integer(dates >= treat_start & dates < other_policy),
    # Seasonality
    summer = as.integer(month %in% 5:9),
    seasonal = sin(2 * pi * (month - 1) / 12),
    seasonal_multiplier = if (is.null(seasonality_index)) {
      exp(0.25 * summer + 0.10 * seasonal)
    } else {
      as.numeric(seasonality_index[month])
    }
  )
}

# =============================================================================
# SPATIAL GRID
# =============================================================================

generate_grid <- function(
  lat_range = c(31, 38),
  lon_range = c(10, 18),
  resolution = 1,
  near_cutoff_km = 200
) {
  # Reference point
  tripoli <- c(lon = 13.19, lat = 32.89)

  # Create grid (cell centroids)
  lat_seq <- seq(lat_range[1], lat_range[2] - resolution, by = resolution)
  lon_seq <- seq(lon_range[1], lon_range[2] - resolution, by = resolution)

  expand.grid(
    lat_center = lat_seq + resolution / 2,
    lon_center = lon_seq + resolution / 2
  ) %>%
    as_tibble() %>%
    mutate(
      unit_id = dplyr::row_number(),
      # Haversine-approximate distance from Tripoli (km)
      dist_tripoli_km = 111.32 * sqrt(
        (lat_center - tripoli["lat"])^2 +
          (cos(lat_center * pi / 180) * (lon_center - tripoli["lon"]))^2
      ),
      dist_scaled = dist_tripoli_km / 100,
      # Binary exposure group for the estimator horse race
      near_libya = as.integer(dist_tripoli_km <= near_cutoff_km)
    )
}

# =============================================================================
# PANEL DATA GENERATION
# =============================================================================

generate_panel_data <- function(
  grid = NULL,
  time_panel = NULL,
  # Treatment effect (log-odds shift in per-attempt mortality when "active")
  delta = 0.6,
  # Baselines
  baseline_attempts = 60,
  baseline_mortality = 0.02,
  # Latent structure (interactive FE)
  n_factors = 2,
  factor_strength = 0.5,
  # NEW: Correlation between loadings and treatment assignment
  # When > 0, near-Libya units have systematically different loadings,
  # violating parallel trends and biasing TWFE
  loading_correlation = 0.5,
  # Measurement issues
  underreporting_rate = 0.0,
  coord_missing_rate = 0.0,
  # Effect timing
  short_lived = FALSE,
  # Grid options
  near_cutoff_km = 200,
  seed = 42
) {
  set.seed(seed)

  if (is.null(grid)) grid <- generate_grid(near_cutoff_km = near_cutoff_km)
  if (is.null(time_panel)) time_panel <- generate_time_panel()

  N <- nrow(grid)
  T_max <- nrow(time_panel)

  # Latent interactive fixed effects
  alpha_i <- rnorm(N, mean = 0, sd = 0.2)
  lambda_t <- rnorm(T_max, mean = 0, sd = 0.15)

  # Factor loadings: random component + systematic shift for near-Libya units
  # This creates correlation between loadings and treatment assignment,
  # which violates parallel trends and biases TWFE
  loadings <- matrix(rnorm(N * n_factors), nrow = N, ncol = n_factors)

  # Add systematic component for treated (near-Libya) units
  # The shift is in the same direction for all factors, creating differential trends
  near_libya_vec <- grid$near_libya
  for (k in 1:n_factors) {
    loadings[near_libya_vec == 1, k] <- loadings[near_libya_vec == 1, k] + loading_correlation
  }

  factors <- matrix(rnorm(T_max * n_factors), nrow = T_max, ncol = n_factors)

  # Make factors have a trend component (not just noise)
  # This ensures the differential loadings translate to differential trends
  for (k in 1:n_factors) {
    trend <- seq(-1, 1, length.out = T_max) * 0.5
    factors[, k] <- factors[, k] + trend * ((-1)^k) # Alternating trend directions
  }

  interactive <- loadings %*% t(factors) # N × T

  # Balanced panel
  panel <- tidyr::crossing(grid, time_panel) %>%
    mutate(
      i = unit_id,
      t = time_id,
      alpha = alpha_i[i],
      lambda = lambda_t[t],
      lf = factor_strength * interactive[cbind(i, t)],

      # Attempts (exposure)
      attempts_mean = baseline_attempts * seasonal_multiplier * exp(
        0.15 * alpha + 0.10 * lambda + 0.10 * lf
      ),
      attempts = rpois(n(), attempts_mean),

      # Treatment indicator used by estimators (binary near×post)
      D = near_libya * post_treat,

      # Effect active: persistent or only in May–June
      effect_active = if (short_lived) near_libya * effect_window else D
    )

  # Mortality probabilities
  panel <- panel %>%
    mutate(
      logit_base = qlogis(baseline_mortality) +
        0.25 * alpha + 0.20 * lambda + 0.30 * lf +
        0.15 * seasonal +
        (-0.015) * dist_scaled,
      p0 = plogis(logit_base),
      p1 = plogis(logit_base + delta * effect_active)
    )

  # Realizations + measurement error + outcomes
  panel <- panel %>%
    mutate(
      deaths_true = rbinom(n(), size = attempts, prob = p1),
      deaths_observed = rbinom(n(), size = deaths_true, prob = 1 - underreporting_rate),
      mortality_rate = ifelse(attempts > 0, deaths_observed / attempts, 0),
      coord_missing = rbinom(n(), 1, coord_missing_rate),
      lat_center = ifelse(coord_missing == 1, NA_real_, lat_center),
      lon_center = ifelse(coord_missing == 1, NA_real_, lon_center),
      Y_any = as.integer(deaths_observed > 0),
      Y_ihs = asinh(deaths_observed),
      EY_any_0 = 1 - (1 - p0)^attempts,
      EY_any_1 = 1 - (1 - p1)^attempts,
      EY_any_obs = 1 - (1 - p1)^attempts
    )

  panel
}

# =============================================================================
# TRUE ESTIMANDS (FOR SIMULATION EVALUATION)
# =============================================================================

compute_true_estimands <- function(panel, outcome_col = "Y_any") {
  if (!("D" %in% names(panel))) {
    panel <- panel %>% mutate(D = near_libya * post_treat)
  }

  treated_post <- panel$D == 1

  if (outcome_col == "Y_any") {
    att <- mean(panel$EY_any_1[treated_post] - panel$EY_any_0[treated_post], na.rm = TRUE)

    tmp <- panel %>%
      mutate(
        group = ifelse(near_libya == 1, "near", "far"),
        period = ifelse(post_treat == 1, "post", "pre")
      ) %>%
      group_by(group, period) %>%
      summarise(m = mean(EY_any_obs, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = period, values_from = m)

    did <- (tmp$post[tmp$group == "near"] - tmp$pre[tmp$group == "near"]) -
      (tmp$post[tmp$group == "far"] - tmp$pre[tmp$group == "far"])

    return(list(att = att, did = did))
  }

  Y <- panel[[outcome_col]]
  tmp <- panel %>%
    mutate(
      group = ifelse(near_libya == 1, "near", "far"),
      period = ifelse(post_treat == 1, "post", "pre")
    ) %>%
    group_by(group, period) %>%
    summarise(m = mean(Y, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = period, values_from = m)

  did <- (tmp$post[tmp$group == "near"] - tmp$pre[tmp$group == "near"]) -
    (tmp$post[tmp$group == "far"] - tmp$pre[tmp$group == "far"])

  list(att = did, did = did)
}

# =============================================================================
# SCENARIOS
# =============================================================================

dgp_scenarios <- function() {
  list(
    # Parallel trends: no latent factors OR loadings uncorrelated with treatment
    "Parallel trends (no latent factors)" = list(
      delta = 0.6, factor_strength = 0.0, loading_correlation = 0.0,
      underreporting_rate = 0.0, short_lived = FALSE
    ),
    # Moderate violation: some latent structure correlated with treatment
    "Moderate latent trends" = list(
      delta = 0.6, factor_strength = 0.5, loading_correlation = 0.5,
      underreporting_rate = 0.0, short_lived = FALSE
    ),
    # Strong violation: strong latent structure highly correlated with treatment
    "Strong latent trends" = list(
      delta = 0.6, factor_strength = 0.9, loading_correlation = 0.8,
      underreporting_rate = 0.0, short_lived = FALSE
    ),
    # Strong violation + measurement error
    "Strong latent + underreporting" = list(
      delta = 0.6, factor_strength = 0.9, loading_correlation = 0.8,
      underreporting_rate = 0.3, short_lived = FALSE
    ),
    # Short-lived effect (transient)
    "Short-lived (May–June only)" = list(
      delta = 0.6, factor_strength = 0.5, loading_correlation = 0.5,
      underreporting_rate = 0.0, short_lived = TRUE
    )
  )
}
