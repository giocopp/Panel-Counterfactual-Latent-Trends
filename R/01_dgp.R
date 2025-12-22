#' ============================================================================
#' Data Generating Process (DGP): Semi-synthetic logistic Interactive FE
#' ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# =============================================================================
# TIME PANEL
# =============================================================================
generate_time_panel <- function(
  start_date = "2015-04-01",
  end_date = "2018-02-01",
  mou_date = "2017-02-01",
  treat_start = "2017-05-01",
  other_policy = "2017-07-01"
) {
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "month")
  mou_date <- as.Date(mou_date)
  treat_start <- as.Date(treat_start)
  other_policy <- as.Date(other_policy)

  tibble(
    date = dates,
    year = year(dates),
    month = month(dates),
    time_id = seq_along(dates),
    post_mou = as.integer(dates >= mou_date),
    post_treat = as.integer(dates >= treat_start),
    effect_window = as.integer(dates >= treat_start & dates < other_policy)
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
  tripoli <- c(lon = 13.19, lat = 32.89)

  lat_seq <- seq(lat_range[1], lat_range[2] - resolution, by = resolution)
  lon_seq <- seq(lon_range[1], lon_range[2] - resolution, by = resolution)

  expand.grid(
    lat_center = lat_seq + resolution / 2,
    lon_center = lon_seq + resolution / 2
  ) %>%
    as_tibble() %>%
    mutate(
      unit_id = dplyr::row_number(),
      unit_key = paste0(floor(lat_center), "_", floor(lon_center)),
      dist_tripoli_km = 111.32 * sqrt(
        (lat_center - tripoli["lat"])^2 +
          (cos(lat_center * pi / 180) * (lon_center - tripoli["lon"]))^2
      ),
      dist_scaled = dist_tripoli_km / 100,
      near_libya = as.integer(dist_tripoli_km <= near_cutoff_km)
    )
}

# =============================================================================
# ANCHOR I/O
# =============================================================================
load_dgp_anchor <- function(path = "data/dgp_anchor.rds") {
  if (!file.exists(path)) stop("Anchor file not found: ", path)
  a <- readRDS(path)
  if (is.null(a$K) || is.null(a$intercept) || is.null(a$alpha) || is.null(a$loadings)) {
    stop("Anchor file is missing required fields (K/intercept/alpha/loadings). Rebuild it.")
  }
  a
}

# =============================================================================
# HELPERS
# =============================================================================
predict_gamma <- function(time_panel, anchor) {
  if (is.null(anchor$gamma_fit) || is.null(anchor$time_center)) {
    stop("Anchor is missing gamma_fit/time_center.")
  }
  newdata <- time_panel %>%
    transmute(
      month = month,
      time_c = time_id - anchor$time_center
    )
  as.numeric(predict(anchor$gamma_fit, newdata = newdata))
}

make_delta_path <- function(time_panel, delta, effect_type = c("constant", "decay", "window"),
                            decay_half_life = 3) {
  effect_type <- match.arg(effect_type)
  post <- time_panel$post_treat

  if (effect_type == "constant") {
    return(delta * post)
  }
  if (effect_type == "window") {
    return(delta * time_panel$effect_window)
  }

  m0 <- min(time_panel$time_id[time_panel$post_treat == 1])
  m <- pmax(0, time_panel$time_id - m0)
  decay <- exp(-log(2) * m / decay_half_life)
  delta * post * decay
}

make_sev_path <- function(time_panel, sev_delta,
                          sev_effect_type = c("decay", "constant", "window"),
                          sev_decay_half_life = 6) {
  sev_effect_type <- match.arg(sev_effect_type)
  post <- time_panel$post_treat

  if (sev_effect_type == "constant") {
    return(sev_delta * post)
  }
  if (sev_effect_type == "window") {
    return(sev_delta * time_panel$effect_window)
  }

  m0 <- min(time_panel$time_id[time_panel$post_treat == 1])
  m <- pmax(0, time_panel$time_id - m0)
  decay <- exp(-log(2) * m / sev_decay_half_life)
  sev_delta * post * decay
}

simulate_factors_ar1 <- function(T_max, phi, sigma, init = NULL) {
  K <- length(phi)
  f <- matrix(0, nrow = T_max, ncol = K)

  for (k in seq_len(K)) {
    ph <- max(min(phi[k], 0.95), -0.95)
    sg <- max(sigma[k], 1e-8)

    if (is.null(init)) {
      f[1, k] <- rnorm(1, 0, sg / sqrt(1 - ph^2))
    } else {
      f[1, k] <- init[k]
    }

    if (T_max >= 2) {
      for (t in 2:T_max) {
        f[t, k] <- ph * f[t - 1, k] + rnorm(1, 0, sg)
      }
    }
  }
  f
}

simulate_u_ar1 <- function(N, T_max, rho, sigma_u) {
  rho <- max(min(rho, 0.95), -0.95)
  sigma_u <- max(sigma_u, 1e-8)

  u <- matrix(0, nrow = N, ncol = T_max)
  u[, 1] <- rnorm(N, 0, sigma_u / sqrt(1 - rho^2))
  if (T_max >= 2) {
    for (t in 2:T_max) {
      u[, t] <- rho * u[, t - 1] + rnorm(N, 0, sigma_u)
    }
  }
  u
}

.expected_asinh_pos <- function(mu, size, n_mc = 5000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- stats::rnbinom(n_mc, size = size, mu = mu)
  x <- pmax(1L, as.integer(x))
  mean(asinh(x))
}

# =============================================================================
# PANEL DATA GENERATION (MAIN ENTRY POINT)
# =============================================================================
generate_panel_data <- function(
  grid = NULL,
  time_panel = NULL,
  anchor = NULL,
  delta = 0.6,
  s = 0.9,
  effect_type = c("constant", "decay", "window"),
  decay_half_life = 3,
  rho = 0.3,
  sigma_u = NULL,
  q = 1.0,
  q_t = NULL,
  auto_level_shift = FALSE,
  target_pre_rate = NULL,
  level_shift = 0,
  mu_pos_cap = Inf,
  sev_delta = 0,
  sev_effect_type = c("decay", "constant", "window"),
  sev_decay_half_life = 6,
  underreporting_rate = NULL,
  short_lived = NULL,
  near_cutoff_km = 200,
  seed = 42
) {
  set.seed(seed)

  if (is.null(grid)) grid <- generate_grid(near_cutoff_km = near_cutoff_km)
  if (is.null(time_panel)) time_panel <- generate_time_panel()
  if (is.null(anchor)) anchor <- load_dgp_anchor("data/dgp_anchor.rds")

  if (!is.null(underreporting_rate)) q <- 1 - underreporting_rate
  if (isTRUE(short_lived)) effect_type <- "window"
  effect_type <- match.arg(effect_type)
  sev_effect_type <- match.arg(sev_effect_type)

  if (is.null(anchor$unit_key)) stop("Anchor is missing unit_key; rebuild anchor.")
  idx <- match(grid$unit_key, anchor$unit_key)
  if (any(is.na(idx))) stop("Grid does not match anchor unit_key. Rebuild anchor for this grid definition.")

  N <- nrow(grid)
  T_max <- nrow(time_panel)
  K <- anchor$K

  intercept <- anchor$intercept
  alpha_i <- as.numeric(anchor$alpha[idx])
  gamma_t <- predict_gamma(time_panel, anchor)

  # Factors + loadings
  loadings <- anchor$loadings[idx, , drop = FALSE]
  phi_f <- anchor$phi_f %||% rep(0.5, K)
  sig_f <- anchor$sigma_f %||% rep(0.2, K)
  f_init <- anchor$f_init %||% rep(0, K)

  factors <- simulate_factors_ar1(T_max, phi = phi_f, sigma = sig_f, init = f_init) # T x K
  lf_mat <- loadings %*% t(factors) # N x T
  lf_mat_scaled <- s * lf_mat # N x T

  # Serial noise
  if (is.null(sigma_u)) sigma_u <- anchor$sigma_u %||% 0.25
  u_mat <- simulate_u_ar1(N, T_max, rho = rho, sigma_u = sigma_u) # N x T

  # Under-reporting vector q_t
  if (is.null(q_t)) {
    q_t <- rep(q, T_max)
  } else {
    if (length(q_t) != T_max) stop("q_t must have length T (nrow(time_panel)).")
  }

  delta_t <- make_delta_path(time_panel, delta = delta, effect_type = effect_type, decay_half_life = decay_half_life)
  sev_t <- make_sev_path(time_panel, sev_delta = sev_delta, sev_effect_type = sev_effect_type, sev_decay_half_life = sev_decay_half_life)

  post_for_D <- if (effect_type == "window") time_panel$effect_window else time_panel$post_treat

  panel <- tidyr::crossing(grid, time_panel) %>%
    mutate(i = unit_id, t = time_id, D = near_libya * post_for_D) %>%
    arrange(unit_id, time_id) %>%
    mutate(
      alpha = alpha_i[unit_id],
      gamma = gamma_t[time_id],
      delta_time = delta_t[time_id],
      treat_effect = near_libya * delta_time,
      sev_time = sev_t[time_id],
      sev_effect = near_libya * sev_time
    )

  # ✅ CRITICAL FIX: correct (unit,time) extraction from N×T matrices
  panel$lf <- as.numeric(lf_mat_scaled[cbind(panel$unit_id, panel$time_id)])
  panel$u <- as.numeric(u_mat[cbind(panel$unit_id, panel$time_id)])

  # auto-calibrate level_shift to match pre-period mean rate
  if (isTRUE(auto_level_shift)) {
    if (is.null(target_pre_rate) || !is.finite(target_pre_rate)) {
      stop("auto_level_shift=TRUE requires target_pre_rate.")
    }
    target_pre_rate <- min(max(target_pre_rate, 1e-6), 1 - 1e-6)

    pre_idx <- which(panel$post_treat == 0)
    base_pre <- intercept + panel$alpha[pre_idx] + panel$gamma[pre_idx] + panel$lf[pre_idx] + panel$u[pre_idx]

    f_root <- function(b) mean(plogis(base_pre + b)) - target_pre_rate
    level_shift <- uniroot(f_root, interval = c(-10, 10))$root
  }

  intercept_adj <- intercept + level_shift

  panel <- panel %>%
    mutate(
      eta0 = intercept_adj + alpha + gamma + lf + u,
      eta1 = eta0 + treat_effect,
      p0 = plogis(eta0),
      p1 = plogis(eta1)
    )

  mu_pos <- anchor$death_mu_pos %||% 3
  size_pos <- anchor$death_size_pos %||% 5
  mu_pos <- min(mu_pos, mu_pos_cap)

  mu_pos_it <- pmin(mu_pos * exp(panel$sev_effect), mu_pos_cap)

  m_asinh0 <- .expected_asinh_pos(mu = mu_pos, size = size_pos, n_mc = 5000, seed = 9001)
  sev_unique <- unique(sev_t)
  m_map <- setNames(
    vapply(sev_unique, function(sd) .expected_asinh_pos(mu = mu_pos * exp(sd), size = size_pos, n_mc = 5000), numeric(1)),
    as.character(sev_unique)
  )
  m_asinh_time <- as.numeric(m_map[as.character(sev_t)])
  m_asinh_it <- m_asinh0 + panel$near_libya * (m_asinh_time[panel$time_id] - m_asinh0)

  panel <- panel %>%
    mutate(
      Y_true = rbinom(n(), size = 1, prob = p1),
      deaths_true = ifelse(
        Y_true == 1,
        pmax(1L, as.integer(stats::rnbinom(n(), size = size_pos, mu = mu_pos_it))),
        0L
      ),
      q_obs = q_t[time_id],
      R_event = rbinom(n(), size = 1, prob = q_obs),
      Y_any = as.integer(Y_true == 1 & R_event == 1),
      deaths_observed = ifelse(Y_any == 1, deaths_true, 0L),
      Y_ihs = asinh(deaths_observed),
      EY_any_0 = p0,
      EY_any_1 = p1,
      EY_any_obs = p1 * q_obs,
      EY_ihs_0 = p0 * q_obs * m_asinh0,
      EY_ihs_1 = p1 * q_obs * m_asinh_it,
      EY_ihs_obs = p1 * q_obs * m_asinh_it,
      level_shift = level_shift
    )

  panel
}

# =============================================================================
# TRUE ESTIMANDS + SCENARIOS (unchanged)
# =============================================================================
compute_true_estimands <- function(
  panel,
  outcome_col = "Y_any",
  treat_type = c("post", "window"),
  trim_after_window = (match.arg(treat_type) == "window")
) {
  treat_type <- match.arg(treat_type)

  if (treat_type == "window" && isTRUE(trim_after_window)) {
    if (!"effect_window" %in% names(panel)) stop("panel must contain effect_window.")
    panel <- panel %>% filter(post_treat == 0 | effect_window == 1)
  }

  post_est <- if (treat_type == "post") panel$post_treat else panel$effect_window
  treated_post <- (panel$near_libya == 1) & (post_est == 1)

  if (outcome_col == "Y_any") {
    att <- mean(panel$EY_any_1[treated_post] - panel$EY_any_0[treated_post], na.rm = TRUE)

    tmp <- panel %>%
      mutate(
        group = ifelse(near_libya == 1, "near", "far"),
        period = ifelse(post_est == 1, "post", "pre")
      ) %>%
      group_by(group, period) %>%
      summarise(m = mean(EY_any_obs, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = period, values_from = m)

    did <- (tmp$post[tmp$group == "near"] - tmp$pre[tmp$group == "near"]) -
      (tmp$post[tmp$group == "far"] - tmp$pre[tmp$group == "far"])

    return(list(att = att, did = did))
  }

  if (outcome_col == "Y_ihs" && all(c("EY_ihs_1", "EY_ihs_0", "EY_ihs_obs") %in% names(panel))) {
    att <- mean(panel$EY_ihs_1[treated_post] - panel$EY_ihs_0[treated_post], na.rm = TRUE)

    tmp <- panel %>%
      mutate(
        group = ifelse(near_libya == 1, "near", "far"),
        period = ifelse(post_est == 1, "post", "pre")
      ) %>%
      group_by(group, period) %>%
      summarise(m = mean(EY_ihs_obs, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = period, values_from = m)

    did <- (tmp$post[tmp$group == "near"] - tmp$pre[tmp$group == "near"]) -
      (tmp$post[tmp$group == "far"] - tmp$pre[tmp$group == "far"])

    return(list(att = att, did = did))
  }

  Y <- panel[[outcome_col]]
  tmp <- panel %>%
    mutate(
      group = ifelse(near_libya == 1, "near", "far"),
      period = ifelse(post_est == 1, "post", "pre")
    ) %>%
    group_by(group, period) %>%
    summarise(m = mean(Y, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = period, values_from = m)

  did <- (tmp$post[tmp$group == "near"] - tmp$pre[tmp$group == "near"]) -
    (tmp$post[tmp$group == "far"] - tmp$pre[tmp$group == "far"])

  list(att = did, did = did)
}

dgp_scenarios <- function(target_pre_rate = NULL, mu_pos_cap = 30) {
  if (is.null(target_pre_rate)) {
    if (file.exists("data/iom_grid_panel.rds")) {
      iom <- readRDS("data/iom_grid_panel.rds")
      target_pre_rate <- mean(iom$Y_any[iom$post_treat == 0])
    } else {
      target_pre_rate <- 0.0664
    }
  }

  list(
    "Baseline (realistic)" = list(
      delta = 0.6, s = 0.9, effect_type = "constant", decay_half_life = 3,
      q = 1.0, auto_level_shift = TRUE, target_pre_rate = target_pre_rate,
      mu_pos_cap = mu_pos_cap
    ),
    "Parallel trends world (s = 0)" = list(
      delta = 0.6, s = 0.0, effect_type = "constant",
      q = 1.0, auto_level_shift = TRUE, target_pre_rate = target_pre_rate,
      mu_pos_cap = mu_pos_cap
    ),
    "Strong latent trends (larger s)" = list(
      delta = 0.6, s = 2, effect_type = "constant",
      q = 1.0, auto_level_shift = TRUE, target_pre_rate = target_pre_rate,
      mu_pos_cap = mu_pos_cap
    ),
    "Short-lived effect (decay)" = list(
      delta = 0.6, s = 0.9, effect_type = "decay", decay_half_life = 2,
      q = 1.0, auto_level_shift = TRUE, target_pre_rate = target_pre_rate,
      mu_pos_cap = mu_pos_cap
    ),
    "Under-reporting (q < 1)" = list(
      delta = 0.6, s = 0.9, effect_type = "constant", q = 0.7,
      auto_level_shift = TRUE, target_pre_rate = target_pre_rate,
      mu_pos_cap = mu_pos_cap
    )
  )
}
