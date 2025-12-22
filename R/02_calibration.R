#' ============================================================================
#' Calibration / Anchoring from IOM Missing Migrants (semi-synthetic DGP)
#'
#' Output: data/dgp_anchor.rds
#'
#' The anchor is estimated using ONLY the pre-period (month < treat_start):
#'   - intercept + alpha_i (unit heterogeneity) on logit scale (smoothed)
#'   - gamma_t (time effects) modeled as factor(month-of-year) + linear trend
#'   - low-rank residual structure via SVD (K factors) -> loadings + AR(1) factors
#'   - residual serial noise AR(1) parameters (rho_u, sigma_u)
#'   - deaths|event distribution (NegBin) for optional count realism
#' ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

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
  iom_cells <- iom_incidents %>%
    mutate(
      lat_cell = floor(lat),
      lon_cell = floor(lon),
      unit_key = paste0(lat_cell, "_", lon_cell)
    )

  grid_key <- grid %>%
    mutate(unit_key = paste0(floor(lat_center), "_", floor(lon_center))) %>%
    select(unit_id, unit_key)

  iom_cells <- iom_cells %>%
    inner_join(grid_key, by = "unit_key")

  cell_month <- iom_cells %>%
    group_by(unit_id, month_date) %>%
    summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop")

  panel <- tidyr::crossing(
    grid %>% select(unit_id),
    time_panel %>% select(date, time_id, month, post_treat, effect_window)
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

# Estimate anchor components from pre-period only
build_dgp_anchor <- function(iom_panel, grid, time_panel,
                             treat_start = "2017-05-01",
                             K = 2,
                             c_smooth = 0.5) {
  treat_start <- as.Date(treat_start)

  # Keys aligned with DGP grid
  grid <- grid %>% arrange(unit_id)
  unit_key <- grid %>%
    mutate(unit_key = paste0(floor(lat_center), "_", floor(lon_center))) %>%
    pull(unit_key)

  pre <- iom_panel %>%
    filter(month_date < treat_start) %>%
    arrange(unit_id, time_id)

  if (nrow(pre) == 0) stop("No pre-period observations after filtering. Check dates/treat_start.")

  # --- Smoothed pre means on probability scale
  # overall
  pbar <- (sum(pre$Y_any) + c_smooth) / (nrow(pre) + 2 * c_smooth)
  pbar <- min(max(pbar, 1e-6), 1 - 1e-6)
  intercept0 <- qlogis(pbar)

  # unit pre means
  unit_stats <- pre %>%
    group_by(unit_id) %>%
    summarise(y = sum(Y_any), Tpre = n(), .groups = "drop") %>%
    mutate(
      p_i = (y + c_smooth) / (Tpre + 2 * c_smooth),
      p_i = pmin(pmax(p_i, 1e-6), 1 - 1e-6),
      alpha_raw = qlogis(p_i) - intercept0
    )

  # time pre means
  time_stats <- pre %>%
    group_by(time_id, month) %>%
    summarise(y = sum(Y_any), N = n(), .groups = "drop") %>%
    mutate(
      p_t = (y + c_smooth) / (N + 2 * c_smooth),
      p_t = pmin(pmax(p_t, 1e-6), 1 - 1e-6),
      gamma_raw = qlogis(p_t) - intercept0
    )

  # --- Recenter alpha/gamma and adjust intercept for identification
  mean_a <- mean(unit_stats$alpha_raw)
  mean_g <- mean(time_stats$gamma_raw)

  intercept <- intercept0 + mean_a + mean_g
  alpha <- unit_stats$alpha_raw - mean_a
  gamma_pre <- time_stats$gamma_raw - mean_g

  # --- Fit gamma model (pre only): factor(month) + linear trend in time_id
  time_center <- mean(time_stats$time_id)
  df_g <- time_stats %>%
    mutate(
      gamma = gamma_pre,
      time_c = time_id - time_center
    )

  gamma_fit <- lm(gamma ~ factor(month) + time_c, data = df_g)

  # --- Build residual matrix on logit scale with smoothing for single Bernoulli
  # cell-month logit proxy: qlogis((Y + c)/(1+2c))
  pre2 <- pre %>%
    left_join(unit_stats %>% select(unit_id) %>% mutate(alpha = alpha), by = "unit_id") %>%
    left_join(time_stats %>% select(time_id, month) %>% mutate(gamma = gamma_pre), by = c("time_id", "month")) %>%
    mutate(
      z = qlogis((Y_any + c_smooth) / (1 + 2 * c_smooth)),
      time_c = time_id - time_center,
      gamma_hat = as.numeric(predict(gamma_fit, newdata = tibble(month = month, time_c = time_c))),
      resid = z - (intercept + alpha + gamma_hat)
    )

  # residual matrix N x Tpre
  Rwide <- pre2 %>%
    select(unit_id, time_id, resid) %>%
    tidyr::pivot_wider(names_from = time_id, values_from = resid) %>%
    arrange(unit_id)

  Rmat <- as.matrix(Rwide %>% select(-unit_id))
  Tpre <- ncol(Rmat)

  K_eff <- min(K, nrow(Rmat), ncol(Rmat))
  if (K_eff < 1) stop("K too small after dimensionality checks.")

  sv <- svd(Rmat, nu = K_eff, nv = K_eff)

  d <- sv$d[seq_len(K_eff)]
  U <- sv$u[, seq_len(K_eff), drop = FALSE]
  V <- sv$v[, seq_len(K_eff), drop = FALSE]

  # Low-rank approx: R ≈ (U sqrt(D)) (V sqrt(D))'
  sqrtD <- diag(sqrt(pmax(d, 0)))
  loadings <- U %*% sqrtD # N x K
  factors_pre <- V %*% sqrtD # Tpre x K

  colnames(loadings) <- paste0("k", seq_len(K_eff))
  colnames(factors_pre) <- paste0("k", seq_len(K_eff))

  # Align loadings to unit_key
  rownames(loadings) <- unit_key[Rwide$unit_id]

  # --- Factor AR(1) params estimated on pre factors
  phi_f <- rep(0.5, K_eff)
  sigma_f <- rep(0.2, K_eff)
  if (Tpre >= 3) {
    for (k in seq_len(K_eff)) {
      f <- factors_pre[, k]
      fit_k <- lm(f[-1] ~ 0 + f[-length(f)])
      ph <- as.numeric(coef(fit_k)[1])
      ph <- max(min(ph, 0.95), -0.95)
      phi_f[k] <- ph
      sigma_f[k] <- sd(residuals(fit_k))
      sigma_f[k] <- max(sigma_f[k], 1e-6)
    }
  }

  # --- Residual after removing rank-K component -> estimate pooled AR(1) noise
  Rhat <- loadings %*% t(factors_pre) # N x Tpre
  E <- Rmat - Rhat

  Elong <- as_tibble(E) %>%
    mutate(unit_id = Rwide$unit_id) %>%
    pivot_longer(cols = -unit_id, names_to = "time_id", values_to = "e") %>%
    mutate(time_id = as.integer(time_id)) %>%
    arrange(unit_id, time_id) %>%
    group_by(unit_id) %>%
    mutate(e_lag = dplyr::lag(e)) %>%
    ungroup() %>%
    filter(!is.na(e_lag))

  rho_u <- 0.3
  sigma_u <- 0.25
  if (nrow(Elong) > 0) {
    fit_u <- lm(e ~ 0 + e_lag, data = Elong)
    rho_u <- as.numeric(coef(fit_u)[1])
    rho_u <- max(min(rho_u, 0.8), -0.8)
    sigma_u <- sd(residuals(fit_u))
    sigma_u <- max(sigma_u, 1e-6)
  }

  # --- Death count model (positive cell-month totals in pre)
  pos <- pre %>% filter(deaths > 0)
  mu_pos <- if (nrow(pos) > 0) mean(pos$deaths) else 3
  var_pos <- if (nrow(pos) > 1) var(pos$deaths) else (mu_pos + 1)

  # NegBin size via method-of-moments: var = mu + mu^2/size  => size = mu^2/(var-mu)
  if (is.na(var_pos) || var_pos <= mu_pos) {
    size_pos <- 1e6
  } else {
    size_pos <- mu_pos^2 / (var_pos - mu_pos)
  }
  size_pos <- max(min(size_pos, 1e6), 0.5)

  # alpha vector aligned to unit_key
  alpha_vec <- alpha
  names(alpha_vec) <- unit_key[unit_stats$unit_id]

  list(
    # meta
    treat_start = treat_start,
    K = K_eff,
    c_smooth = c_smooth,

    # identification + FE
    intercept = intercept,
    alpha = alpha_vec,

    # gamma model
    gamma_fit = gamma_fit,
    time_center = time_center,

    # low rank
    unit_key = unit_key,
    loadings = loadings,
    phi_f = phi_f,
    sigma_f = sigma_f,
    f_init = as.numeric(factors_pre[1, , drop = TRUE]),

    # serial noise
    rho_u = rho_u,
    sigma_u = sigma_u,

    # deaths model
    death_mu_pos = mu_pos,
    death_size_pos = size_pos
  )
}

# Main entry point used by run scripts
calibrate_from_iom <- function(
  iom_path = NULL,
  lat_range = c(31, 38),
  lon_range = c(10, 18),
  start_date = "2015-04-01",
  end_date = "2018-02-01",
  treat_start = "2017-05-01",
  near_cutoff_km = 200,
  K = 2,
  c_smooth = 0.5
) {
  if (is.null(iom_path)) iom_path <- find_iom_file()
  if (is.null(iom_path)) {
    message("No IOM incident file found; cannot build anchor.")
    return(NULL)
  }

  # Use the same grid/time panel definitions as the DGP
  grid <- generate_grid(lat_range = lat_range, lon_range = lon_range, near_cutoff_km = near_cutoff_km)
  time_panel <- generate_time_panel(
    start_date = start_date,
    end_date = end_date,
    treat_start = treat_start
  )

  iom <- read_iom_incidents(iom_path)
  iom_panel <- build_iom_grid_panel(iom, grid, time_panel)

  anchor <- build_dgp_anchor(
    iom_panel = iom_panel,
    grid = grid,
    time_panel = time_panel,
    treat_start = treat_start,
    K = K,
    c_smooth = c_smooth
  )

  dir.create("data", showWarnings = FALSE, recursive = TRUE)
  saveRDS(anchor, "data/dgp_anchor.rds")
  saveRDS(iom_panel, "data/iom_grid_panel.rds")

  message("Saved anchor to data/dgp_anchor.rds")
  message("Saved mapped IOM panel to data/iom_grid_panel.rds")

  list(
    iom_path = iom_path,
    grid = grid,
    time_panel = time_panel,
    iom_panel = iom_panel,
    anchor = anchor
  )
}
