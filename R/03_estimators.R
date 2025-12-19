#' ============================================================================
#' Estimators: TWFE, Matrix Completion, Synthetic DiD, TROP
#'
#' Target estimand:
#'   D_it = near_libya_i × post_treat_t
#'   ATT on a chosen outcome Y (default: Y_any)
#'
#' Notes:
#' - Y_any is binary; we treat it as continuous for comparability across methods
#'   (linear probability / least-squares imputation).
#' ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
  library(Matrix)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

prepare_estimation_data <- function(panel, outcome_col = "Y_any") {
  stopifnot(outcome_col %in% names(panel))
  panel %>%
    mutate(
      D = near_libya * post_treat,
      Y = .data[[outcome_col]]
    )
}

# =============================================================================
# TWFE
# =============================================================================

estimate_twfe <- function(panel, outcome_col = "Y_any", cluster = "unit_id") {
  data <- prepare_estimation_data(panel, outcome_col = outcome_col)
  fit <- feols(Y ~ D | unit_id + time_id, data = data, cluster = cluster)

  ct <- as.data.frame(fixest::coeftable(fit))
  est <- ct["D", "Estimate"]
  se  <- ct["D", "Std. Error"]
  pval <- ct["D", "Pr(>|t|)"] %||% ct["D", "Pr(>|z|)"]

  list(
    method = "TWFE",
    estimand = outcome_col,
    estimate = as.numeric(est),
    se = as.numeric(se),
    ci_lower = as.numeric(est - 1.96 * se),
    ci_upper = as.numeric(est + 1.96 * se),
    p_value = as.numeric(pval),
    n_obs = fit$nobs,
    fit = fit
  )
}

# Optional, for exposition: continuous distance×post TWFE
estimate_twfe_distance_post <- function(panel, outcome_col = "Y_any", cluster = "unit_id") {
  stopifnot("dist_tripoli_km" %in% names(panel))
  data <- prepare_estimation_data(panel, outcome_col = outcome_col)

  fit <- feols(
    Y ~ dist_tripoli_km:post_treat + dist_tripoli_km + post_treat | unit_id + time_id,
    data = data, cluster = cluster
  )

  list(method = "TWFE (distance×post)", fit = fit)
}

# =============================================================================
# MATRIX COMPLETION (simple SVD completion)
# =============================================================================

estimate_mc <- function(panel, outcome_col = "Y_any", lambda = 1, max_rank = 3, n_boot = 50) {
  data <- prepare_estimation_data(panel, outcome_col = outcome_col)

  Y_wide <- data %>%
    select(unit_id, time_id, Y) %>%
    pivot_wider(names_from = time_id, values_from = Y) %>%
    arrange(unit_id) %>%
    select(-unit_id) %>%
    as.matrix()

  D_wide <- data %>%
    select(unit_id, time_id, D) %>%
    pivot_wider(names_from = time_id, values_from = D) %>%
    arrange(unit_id) %>%
    select(-unit_id) %>%
    as.matrix()

  N <- nrow(Y_wide)

  Y_masked <- Y_wide
  Y_masked[D_wide == 1] <- NA

  row_means <- rowMeans(Y_masked, na.rm = TRUE)
  row_means[is.nan(row_means)] <- 0
  for (i in 1:N) {
    Y_masked[i, is.na(Y_masked[i, ])] <- row_means[i]
  }

  svd_fit <- svd(Y_masked)
  d_thresh <- pmax(svd_fit$d - lambda, 0)
  k <- min(max_rank, sum(d_thresh > 0))
  if (k == 0) k <- 1

  Y_completed <- svd_fit$u[, 1:k, drop = FALSE] %*%
    diag(d_thresh[1:k], nrow = k) %*%
    t(svd_fit$v[, 1:k, drop = FALSE])

  treated_entries <- D_wide == 1
  att <- mean(Y_wide[treated_entries] - Y_completed[treated_entries], na.rm = TRUE)

  att_boot <- rep(NA_real_, n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(N), N, replace = TRUE)
    Y_b <- Y_wide[idx, , drop = FALSE]
    D_b <- D_wide[idx, , drop = FALSE]
    Y_cf_b <- Y_completed[idx, , drop = FALSE]
    treated_b <- D_b == 1
    att_boot[b] <- mean(Y_b[treated_b] - Y_cf_b[treated_b], na.rm = TRUE)
  }
  se <- sd(att_boot, na.rm = TRUE)

  list(
    method = "Matrix Completion",
    estimand = outcome_col,
    estimate = as.numeric(att),
    se = as.numeric(se),
    ci_lower = as.numeric(att - 1.96 * se),
    ci_upper = as.numeric(att + 1.96 * se),
    p_value = as.numeric(2 * pnorm(-abs(att / se))),
    n_obs = sum(!is.na(Y_wide)),
    Y_completed = Y_completed
  )
}

# =============================================================================
# SYNTHETIC DiD (simplified)
# =============================================================================

estimate_sdid <- function(panel, outcome_col = "Y_any") {
  data <- prepare_estimation_data(panel, outcome_col = outcome_col)

  Y_wide <- data %>%
    select(unit_id, time_id, Y) %>%
    pivot_wider(names_from = time_id, values_from = Y) %>%
    arrange(unit_id) %>%
    select(-unit_id) %>%
    as.matrix()

  unit_treated <- data %>%
    group_by(unit_id) %>%
    summarise(treated = max(near_libya), .groups = "drop") %>%
    arrange(unit_id) %>%
    pull(treated)

  time_info <- data %>%
    distinct(time_id, post_treat) %>%
    arrange(time_id)

  pre_periods  <- time_info$time_id[time_info$post_treat == 0]
  post_periods <- time_info$time_id[time_info$post_treat == 1]

  idx_treat <- which(unit_treated == 1)
  idx_ctrl  <- which(unit_treated == 0)

  Y_treat <- Y_wide[idx_treat, , drop = FALSE]
  Y_ctrl  <- Y_wide[idx_ctrl,  , drop = FALSE]

  Y_ctrl_pre <- Y_ctrl[, pre_periods, drop = FALSE]
  Y_treat_pre_mean <- colMeans(Y_treat[, pre_periods, drop = FALSE], na.rm = TRUE)

  Y_ctrl_pre_filled <- Y_ctrl_pre
  for (j in seq_len(ncol(Y_ctrl_pre_filled))) {
    colm <- mean(Y_ctrl_pre_filled[, j], na.rm = TRUE)
    Y_ctrl_pre_filled[is.na(Y_ctrl_pre_filled[, j]), j] <- colm
  }

  distances <- rowSums((Y_ctrl_pre_filled -
                          matrix(Y_treat_pre_mean, nrow(Y_ctrl_pre_filled), length(pre_periods), byrow = TRUE))^2,
                       na.rm = TRUE)
  omega <- exp(-distances / (2 * var(distances) + 1e-6))
  omega <- omega / sum(omega)

  lambda_t <- rep(1 / length(pre_periods), length(pre_periods))

  Y_treat_post <- mean(Y_treat[, post_periods], na.rm = TRUE)
  Y_treat_pre_w <- sum(lambda_t * colMeans(Y_treat[, pre_periods], na.rm = TRUE))

  Y_ctrl_post_w <- sum(omega * rowMeans(Y_ctrl[, post_periods, drop = FALSE], na.rm = TRUE))
  Y_ctrl_pre_w  <- sum(omega * (Y_ctrl_pre_filled %*% lambda_t))

  att <- (Y_treat_post - Y_treat_pre_w) - (Y_ctrl_post_w - Y_ctrl_pre_w)

  n_treat <- nrow(Y_treat)
  att_jack <- rep(NA_real_, n_treat)
  for (i in seq_len(n_treat)) {
    Y_t_loo <- Y_treat[-i, , drop = FALSE]
    Y_t_post_loo <- mean(Y_t_loo[, post_periods], na.rm = TRUE)
    Y_t_pre_loo  <- sum(lambda_t * colMeans(Y_t_loo[, pre_periods], na.rm = TRUE))
    att_jack[i] <- (Y_t_post_loo - Y_t_pre_loo) - (Y_ctrl_post_w - Y_ctrl_pre_w)
  }
  se <- sqrt((n_treat - 1) / n_treat * sum((att_jack - mean(att_jack, na.rm = TRUE))^2, na.rm = TRUE))

  list(
    method = "Synthetic DiD",
    estimand = outcome_col,
    estimate = as.numeric(att),
    se = as.numeric(se),
    ci_lower = as.numeric(att - 1.96 * se),
    ci_upper = as.numeric(att + 1.96 * se),
    p_value = as.numeric(2 * pnorm(-abs(att / se))),
    n_obs = sum(!is.na(Y_wide))
  )
}

# =============================================================================
# TROP (simplified hybrid)
# =============================================================================

estimate_trop <- function(
    panel,
    outcome_col = "Y_any",
    max_rank = 2,
    lambda_L = 1,
    n_boot = 50
) {
  data <- prepare_estimation_data(panel, outcome_col = outcome_col)

  Y_wide <- data %>%
    select(unit_id, time_id, Y) %>%
    pivot_wider(names_from = time_id, values_from = Y) %>%
    arrange(unit_id) %>%
    select(-unit_id) %>%
    as.matrix()

  D_wide <- data %>%
    select(unit_id, time_id, D) %>%
    pivot_wider(names_from = time_id, values_from = D) %>%
    arrange(unit_id) %>%
    select(-unit_id) %>%
    as.matrix()

  unit_treated <- data %>%
    group_by(unit_id) %>%
    summarise(treated = max(near_libya), .groups = "drop") %>%
    arrange(unit_id) %>%
    pull(treated)

  time_info <- data %>%
    distinct(time_id, post_treat) %>%
    arrange(time_id)

  N <- nrow(Y_wide)
  T_obs <- ncol(Y_wide)
  pre_idx  <- which(time_info$post_treat == 0)
  post_idx <- which(time_info$post_treat == 1)
  treat_idx <- which(unit_treated == 1)
  ctrl_idx  <- which(unit_treated == 0)

  # Step 1: TWFE + low-rank on controls (mask treated entries)
  Y_ctrl <- Y_wide
  Y_ctrl[D_wide == 1] <- NA

  row_means <- rowMeans(Y_ctrl, na.rm = TRUE); row_means[is.nan(row_means)] <- 0
  col_means <- colMeans(Y_ctrl, na.rm = TRUE); col_means[is.nan(col_means)] <- 0
  grand_mean <- mean(Y_ctrl, na.rm = TRUE)

  Y_twfe <- outer(row_means, rep(1, T_obs)) + outer(rep(1, N), col_means) - grand_mean
  Y_resid <- Y_ctrl - Y_twfe
  Y_resid[is.na(Y_resid)] <- 0

  svd_fit <- svd(Y_resid)
  d_thresh <- pmax(svd_fit$d - lambda_L, 0)
  k <- min(max_rank, sum(d_thresh > 0))
  if (k == 0) k <- 1

  L_hat <- svd_fit$u[, 1:k, drop = FALSE] %*%
    diag(d_thresh[1:k], nrow = k) %*%
    t(svd_fit$v[, 1:k, drop = FALSE])
  Y_hat <- Y_twfe + L_hat

  # TROP ATT: compare treated observed to control-predicted trajectory
  Y_treat_post <- mean(Y_wide[treat_idx, post_idx, drop = FALSE], na.rm = TRUE)
  Y_treat_pre  <- mean(Y_wide[treat_idx, pre_idx,  drop = FALSE], na.rm = TRUE)

  Y_ctrl_post_hat <- mean(Y_hat[ctrl_idx, post_idx, drop = FALSE], na.rm = TRUE)
  Y_ctrl_pre_hat  <- mean(Y_hat[ctrl_idx, pre_idx,  drop = FALSE], na.rm = TRUE)

  att <- (Y_treat_post - Y_treat_pre) - (Y_ctrl_post_hat - Y_ctrl_pre_hat)

  # Bootstrap SE by resampling units (conservative proxy)
  att_boot <- rep(NA_real_, n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(N), N, replace = TRUE)
    Yb <- Y_wide[idx, , drop = FALSE]
    treated_b <- unit_treated[idx]
    tb <- which(treated_b == 1)
    cb <- which(treated_b == 0)
    if (length(tb) < 2 || length(cb) < 2) next

    Y_treat_post_b <- mean(Yb[tb, post_idx, drop = FALSE], na.rm = TRUE)
    Y_treat_pre_b  <- mean(Yb[tb, pre_idx,  drop = FALSE], na.rm = TRUE)
    Y_ctrl_post_b  <- mean(Y_hat[idx, , drop = FALSE][cb, post_idx, drop = FALSE], na.rm = TRUE)
    Y_ctrl_pre_b   <- mean(Y_hat[idx, , drop = FALSE][cb, pre_idx,  drop = FALSE], na.rm = TRUE)
    att_boot[b] <- (Y_treat_post_b - Y_treat_pre_b) - (Y_ctrl_post_b - Y_ctrl_pre_b)
  }
  se <- sd(att_boot, na.rm = TRUE)

  list(
    method = "TROP",
    estimand = outcome_col,
    estimate = as.numeric(att),
    se = as.numeric(se),
    ci_lower = as.numeric(att - 1.96 * se),
    ci_upper = as.numeric(att + 1.96 * se),
    p_value = as.numeric(2 * pnorm(-abs(att / se))),
    n_obs = sum(!is.na(Y_wide))
  )
}

# =============================================================================
# RUN ALL
# =============================================================================

run_all_estimators <- function(panel, outcome_col = "Y_any", true_att = NULL) {
  results <- list(
    twfe = tryCatch(estimate_twfe(panel, outcome_col = outcome_col),
                    error = function(e) list(method = "TWFE", estimate = NA_real_, se = NA_real_, p_value = NA_real_,
                                             ci_lower = NA_real_, ci_upper = NA_real_)),
    mc = tryCatch(estimate_mc(panel, outcome_col = outcome_col),
                  error = function(e) list(method = "Matrix Completion", estimate = NA_real_, se = NA_real_, p_value = NA_real_,
                                           ci_lower = NA_real_, ci_upper = NA_real_)),
    sdid = tryCatch(estimate_sdid(panel, outcome_col = outcome_col),
                    error = function(e) list(method = "Synthetic DiD", estimate = NA_real_, se = NA_real_, p_value = NA_real_,
                                             ci_lower = NA_real_, ci_upper = NA_real_)),
    trop = tryCatch(estimate_trop(panel, outcome_col = outcome_col),
                    error = function(e) list(method = "TROP", estimate = NA_real_, se = NA_real_, p_value = NA_real_,
                                             ci_lower = NA_real_, ci_upper = NA_real_))
  )

  results_df <- purrr::map_dfr(results, function(r) {
    tibble(
      method = r$method,
      estimand = outcome_col,
      estimate = r$estimate,
      se = r$se,
      p_value = r$p_value,
      ci_lower = r$ci_lower %||% (r$estimate - 1.96 * r$se),
      ci_upper = r$ci_upper %||% (r$estimate + 1.96 * r$se),
      n_obs = r$n_obs %||% NA_integer_
    )
  })

  if (!is.null(true_att)) {
    results_df <- results_df %>%
      mutate(
        true_att = true_att,
        bias = estimate - true_att,
        covers = ci_lower <= true_att & ci_upper >= true_att,
        rejects = p_value < 0.05
      )
  }
  results_df
}
