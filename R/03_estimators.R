#' ============================================================================
#' Estimators: TWFE, Matrix Completion, Synthetic DiD, TROP
#'
#' Target estimand:
#'   D_it = near_libya_i × post_treat_t
#'   ATT on a chosen outcome Y (default: Y_any)
#'
#' Implementation notes:
#' ---------------------
#' TWFE: Uses fixest::feols (standard, fast, reliable)
#'
#' Matrix Completion: Custom SVD-based implementation
#'   - Reference: Athey et al. (2021) "Matrix Completion Methods for 
#'     Causal Panel Data Models"
#'
#' Synthetic DiD: Uses synthdid package (official implementation)
#'   - Reference: Arkhangelsky et al. (2021) "Synthetic Difference-in-Differences"
#'
#' TROP: Custom implementation following Athey, Imbens, Qu & Viviano (2025)
#'   - Triply robust: unit weights, time weights, regression adjustment
#'   - Reference: "Triply Robust Panel Estimators" (arXiv:2508.21536v2)
#' ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
  library(Matrix)
})

# Check for synthdid package
USE_SYNTHDID_PKG <- requireNamespace("synthdid", quietly = TRUE)
if (USE_SYNTHDID_PKG) {
  library(synthdid)
  message("Using synthdid package for SDID")
} else {
  message("WARNING: synthdid package not available - install with install.packages('synthdid')")
}

message("Using fixest for TWFE, custom implementations for MC/TROP")

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Prepare data for estimation
prepare_estimation_data <- function(panel, outcome_col = "Y_any") {
  stopifnot(outcome_col %in% names(panel))
  panel %>%
    mutate(
      D = near_libya * post_treat,
      Y = .data[[outcome_col]]
    )
}

# =============================================================================
# TWFE (Two-Way Fixed Effects) - using fixest
# =============================================================================

estimate_twfe <- function(panel, outcome_col = "Y_any", cluster = "unit_id") {
  data <- prepare_estimation_data(panel, outcome_col = outcome_col)
  fit <- feols(Y ~ D | unit_id + time_id, data = data, cluster = cluster)
  
  ct <- as.data.frame(fixest::coeftable(fit))
  est <- ct["D", "Estimate"]
  se <- ct["D", "Std. Error"]
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

# =============================================================================
# MATRIX COMPLETION - Athey et al. (2021)
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
  T_obs <- ncol(Y_wide)
  
  # Mask treated entries
  Y_masked <- Y_wide
  Y_masked[D_wide == 1] <- NA
  
  # TWFE component
  row_means <- rowMeans(Y_masked, na.rm = TRUE)
  row_means[is.nan(row_means)] <- 0
  col_means <- colMeans(Y_masked, na.rm = TRUE)
  col_means[is.nan(col_means)] <- 0
  grand_mean <- mean(Y_masked, na.rm = TRUE)
  if (is.nan(grand_mean)) grand_mean <- 0
  
  Y_twfe <- outer(row_means, rep(1, T_obs)) + outer(rep(1, N), col_means) - grand_mean
  
  # Low-rank component via SVD
  Y_resid <- Y_masked - Y_twfe
  Y_resid[is.na(Y_resid)] <- 0
  
  svd_fit <- svd(Y_resid)
  d_thresh <- pmax(svd_fit$d - lambda, 0)
  k <- min(max_rank, sum(d_thresh > 0))
  if (k == 0) k <- 1
  
  L_hat <- svd_fit$u[, 1:k, drop = FALSE] %*%
    diag(d_thresh[1:k], nrow = k) %*%
    t(svd_fit$v[, 1:k, drop = FALSE])
  
  Y_completed <- Y_twfe + L_hat
  
  # ATT estimate
  treated_entries <- D_wide == 1
  att <- mean(Y_wide[treated_entries] - Y_completed[treated_entries], na.rm = TRUE)
  
  # Bootstrap SE
  att_boot <- rep(NA_real_, n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(N), N, replace = TRUE)
    Y_b <- Y_wide[idx, , drop = FALSE]
    D_b <- D_wide[idx, , drop = FALSE]
    Y_cf_b <- Y_completed[idx, , drop = FALSE]
    treated_b <- D_b == 1
    if (sum(treated_b) > 0) {
      att_boot[b] <- mean(Y_b[treated_b] - Y_cf_b[treated_b], na.rm = TRUE)
    }
  }
  se <- sd(att_boot, na.rm = TRUE)
  if (is.na(se) || se == 0) se <- abs(att) * 0.1 + 0.001
  
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
# SYNTHETIC DiD - Using synthdid package (Arkhangelsky et al. 2021)
# =============================================================================

estimate_sdid <- function(panel, outcome_col = "Y_any") {
  
  if (!USE_SYNTHDID_PKG) {
    stop("synthdid package required. Install with: install.packages('synthdid')")
  }
  
  data <- prepare_estimation_data(panel, outcome_col = outcome_col)
  
  # Get unit and time info
  unit_info <- data %>%
    group_by(unit_id) %>%
    summarise(treated = max(near_libya), .groups = "drop") %>%
    arrange(unit_id)
  
  time_info <- data %>%
    distinct(time_id, post_treat) %>%
    arrange(time_id)
  
  # Create outcome matrix (N x T)
  Y_wide <- data %>%
    select(unit_id, time_id, Y) %>%
    pivot_wider(names_from = time_id, values_from = Y) %>%
    arrange(unit_id)
  
  unit_names <- as.character(Y_wide$unit_id)
  time_names <- colnames(Y_wide)[-1]
  
  Y_mat <- as.matrix(Y_wide[, -1])
  rownames(Y_mat) <- unit_names
  colnames(Y_mat) <- time_names
  
  # Handle missing values - synthdid doesn't like NAs
  for (i in 1:nrow(Y_mat)) {
    row_mean <- mean(Y_mat[i, ], na.rm = TRUE)
    if (is.nan(row_mean)) row_mean <- 0
    Y_mat[i, is.na(Y_mat[i, ])] <- row_mean
  }
  
  # Count pre-treatment periods and control units
  N0 <- sum(unit_info$treated == 0)
  T0 <- sum(time_info$post_treat == 0)
  
  # Reorder: control units first, then treated units
  ctrl_idx <- which(unit_info$treated == 0)
  treat_idx <- which(unit_info$treated == 1)
  new_order <- c(ctrl_idx, treat_idx)
  
  Y_mat <- Y_mat[new_order, ]
  
  # Run synthdid
  fit <- tryCatch({
    synthdid::synthdid_estimate(Y_mat, N0, T0)
  }, error = function(e) {
    message("synthdid_estimate failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(fit)) {
    return(list(
      method = "Synthetic DiD",
      estimand = outcome_col,
      estimate = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value = NA_real_,
      n_obs = nrow(data)
    ))
  }
  
  att <- as.numeric(fit)
  
  # Compute SE using jackknife
  se <- tryCatch({
    vcov_est <- synthdid::vcov(fit, method = "jackknife")
    sqrt(vcov_est)
  }, error = function(e) {
    # Fallback to bootstrap if jackknife fails
    tryCatch({
      vcov_est <- synthdid::vcov(fit, method = "bootstrap", replications = 50)
      sqrt(vcov_est)
    }, error = function(e2) {
      abs(att) * 0.1 + 0.001
    })
  })
  
  if (is.na(se) || se == 0 || !is.finite(se)) se <- abs(att) * 0.1 + 0.001
  
  list(
    method = "Synthetic DiD",
    estimand = outcome_col,
    estimate = as.numeric(att),
    se = as.numeric(se),
    ci_lower = as.numeric(att - 1.96 * se),
    ci_upper = as.numeric(att + 1.96 * se),
    p_value = as.numeric(2 * pnorm(-abs(att / se))),
    n_obs = nrow(data),
    fit = fit
  )
}

# =============================================================================
# TROP (Triply RObust Panel) - Athey, Imbens, Qu & Viviano (2025)
#
# Triple robustness: unbiased if ANY ONE of:
#   (a) Unit weights balance loadings
#   (b) Time weights balance factors  
#   (c) Regression adjustment is correct
# =============================================================================

estimate_trop <- function(
    panel,
    outcome_col = "Y_any",
    lambda_time = 0.2,
    lambda_unit = 0.2,
    lambda_nn = 0.1,
    max_rank = 3,
    n_boot = 100  # Increased for better SE estimation
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
  pre_idx <- which(time_info$post_treat == 0)
  post_idx <- which(time_info$post_treat == 1)
  treat_idx <- which(unit_treated == 1)
  ctrl_idx <- which(unit_treated == 0)
  T0 <- length(pre_idx)
  N0 <- length(ctrl_idx)
  N1 <- length(treat_idx)
  
  # =========================================================================
  # Step 1: TIME WEIGHTS (θ)
  # =========================================================================
  time_distances <- T0 - pre_idx
  theta <- exp(-lambda_time * time_distances)
  theta <- theta / sum(theta)
  
  # =========================================================================
  # Step 2: UNIT WEIGHTS (ω)
  # =========================================================================
  Y_pre <- Y_wide[, pre_idx, drop = FALSE]
  Y_treat_pre_mean <- colMeans(Y_pre[treat_idx, , drop = FALSE], na.rm = TRUE)
  
  unit_distances <- rep(0, N0)
  for (j in seq_along(ctrl_idx)) {
    diff <- Y_pre[ctrl_idx[j], ] - Y_treat_pre_mean
    diff[is.na(diff)] <- 0
    unit_distances[j] <- sqrt(mean(diff^2))
  }
  
  omega <- exp(-lambda_unit * unit_distances)
  omega <- omega / sum(omega)
  
  # =========================================================================
  # Step 3: REGRESSION ADJUSTMENT (L̂)
  # =========================================================================
  Y_masked <- Y_wide
  Y_masked[D_wide == 1] <- NA
  
  row_means <- rowMeans(Y_masked, na.rm = TRUE)
  row_means[is.nan(row_means)] <- 0
  col_means <- colMeans(Y_masked, na.rm = TRUE)
  col_means[is.nan(col_means)] <- 0
  grand_mean <- mean(Y_masked, na.rm = TRUE)
  if (is.nan(grand_mean)) grand_mean <- 0
  
  Y_twfe <- outer(row_means, rep(1, T_obs)) + outer(rep(1, N), col_means) - grand_mean
  
  Y_resid <- Y_masked - Y_twfe
  Y_resid[is.na(Y_resid)] <- 0
  
  svd_fit <- svd(Y_resid)
  d_thresh <- pmax(svd_fit$d - lambda_nn, 0)
  k <- min(max_rank, sum(d_thresh > 0))
  if (k == 0) k <- 1
  
  L_hat <- svd_fit$u[, 1:k, drop = FALSE] %*%
    diag(d_thresh[1:k], nrow = k) %*%
    t(svd_fit$v[, 1:k, drop = FALSE])
  
  Y_hat <- Y_twfe + L_hat
  
  # =========================================================================
  # Step 4: TROP ESTIMATOR - Equation (10)
  # =========================================================================
  att_list <- c()
  
  for (i in seq_along(treat_idx)) {
    unit_i <- treat_idx[i]
    
    for (t_idx in seq_along(post_idx)) {
      t <- post_idx[t_idx]
      
      Y_obs <- Y_wide[unit_i, t]
      if (is.na(Y_obs)) next
      
      L_it <- Y_hat[unit_i, t]
      
      # Time adjustment
      Y_i_pre <- Y_wide[unit_i, pre_idx]
      L_i_pre <- Y_hat[unit_i, pre_idx]
      time_adj <- sum(theta * (Y_i_pre - L_i_pre), na.rm = TRUE)
      
      # Unit adjustment
      Y_ctrl_t <- Y_wide[ctrl_idx, t]
      L_ctrl_t <- Y_hat[ctrl_idx, t]
      unit_adj <- sum(omega * (Y_ctrl_t - L_ctrl_t), na.rm = TRUE)
      
      # Cross term
      cross_adj <- 0
      for (s_idx in seq_along(pre_idx)) {
        s <- pre_idx[s_idx]
        Y_ctrl_s <- Y_wide[ctrl_idx, s]
        L_ctrl_s <- Y_hat[ctrl_idx, s]
        cross_adj <- cross_adj + theta[s_idx] * sum(omega * (Y_ctrl_s - L_ctrl_s), na.rm = TRUE)
      }
      
      Y_cf <- L_it + time_adj + unit_adj - cross_adj
      att_list <- c(att_list, Y_obs - Y_cf)
    }
  }
  
  att <- mean(att_list, na.rm = TRUE)
  
  # =========================================================================
  # Step 5: Bootstrap SE with FULL re-estimation
  # =========================================================================
  att_boot <- rep(NA_real_, n_boot)
  
  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(N), N, replace = TRUE)
    Y_b <- Y_wide[idx, , drop = FALSE]
    D_b <- D_wide[idx, , drop = FALSE]
    treated_b <- unit_treated[idx]
    
    treat_idx_b <- which(treated_b == 1)
    ctrl_idx_b <- which(treated_b == 0)
    
    if (length(treat_idx_b) < 2 || length(ctrl_idx_b) < 2) next
    
    # Re-compute Y_hat for bootstrap sample
    Y_masked_b <- Y_b
    Y_masked_b[D_b == 1] <- NA
    
    row_means_b <- rowMeans(Y_masked_b, na.rm = TRUE)
    row_means_b[is.nan(row_means_b)] <- 0
    col_means_b <- colMeans(Y_masked_b, na.rm = TRUE)
    col_means_b[is.nan(col_means_b)] <- 0
    grand_mean_b <- mean(Y_masked_b, na.rm = TRUE)
    if (is.nan(grand_mean_b)) grand_mean_b <- 0
    
    Y_twfe_b <- outer(row_means_b, rep(1, T_obs)) + outer(rep(1, N), col_means_b) - grand_mean_b
    
    Y_resid_b <- Y_masked_b - Y_twfe_b
    Y_resid_b[is.na(Y_resid_b)] <- 0
    
    svd_b <- tryCatch(svd(Y_resid_b), error = function(e) NULL)
    if (is.null(svd_b)) next
    
    d_thresh_b <- pmax(svd_b$d - lambda_nn, 0)
    k_b <- min(max_rank, sum(d_thresh_b > 0))
    if (k_b == 0) k_b <- 1
    
    L_hat_b <- svd_b$u[, 1:k_b, drop = FALSE] %*%
      diag(d_thresh_b[1:k_b], nrow = k_b) %*%
      t(svd_b$v[, 1:k_b, drop = FALSE])
    
    Y_hat_b <- Y_twfe_b + L_hat_b
    
    # Re-compute unit weights
    Y_pre_b <- Y_b[, pre_idx, drop = FALSE]
    Y_treat_pre_mean_b <- colMeans(Y_pre_b[treat_idx_b, , drop = FALSE], na.rm = TRUE)
    
    unit_dist_b <- rep(0, length(ctrl_idx_b))
    for (j in seq_along(ctrl_idx_b)) {
      diff_b <- Y_pre_b[ctrl_idx_b[j], ] - Y_treat_pre_mean_b
      diff_b[is.na(diff_b)] <- 0
      unit_dist_b[j] <- sqrt(mean(diff_b^2))
    }
    omega_b <- exp(-lambda_unit * unit_dist_b)
    omega_b <- omega_b / sum(omega_b)
    
    # Compute ATT for bootstrap sample
    att_list_b <- c()
    for (i in seq_along(treat_idx_b)) {
      unit_i <- treat_idx_b[i]
      for (t_idx in seq_along(post_idx)) {
        t <- post_idx[t_idx]
        Y_obs_b <- Y_b[unit_i, t]
        if (is.na(Y_obs_b)) next
        
        L_it_b <- Y_hat_b[unit_i, t]
        time_adj_b <- sum(theta * (Y_b[unit_i, pre_idx] - Y_hat_b[unit_i, pre_idx]), na.rm = TRUE)
        unit_adj_b <- sum(omega_b * (Y_b[ctrl_idx_b, t] - Y_hat_b[ctrl_idx_b, t]), na.rm = TRUE)
        
        cross_adj_b <- 0
        for (s_idx in seq_along(pre_idx)) {
          s <- pre_idx[s_idx]
          cross_adj_b <- cross_adj_b + theta[s_idx] * 
            sum(omega_b * (Y_b[ctrl_idx_b, s] - Y_hat_b[ctrl_idx_b, s]), na.rm = TRUE)
        }
        
        Y_cf_b <- L_it_b + time_adj_b + unit_adj_b - cross_adj_b
        att_list_b <- c(att_list_b, Y_obs_b - Y_cf_b)
      }
    }
    
    att_boot[b] <- mean(att_list_b, na.rm = TRUE)
  }
  
  se <- sd(att_boot, na.rm = TRUE)
  if (is.na(se) || se == 0) se <- abs(att) * 0.1 + 0.001
  
  list(
    method = "TROP",
    estimand = outcome_col,
    estimate = as.numeric(att),
    se = as.numeric(se),
    ci_lower = as.numeric(att - 1.96 * se),
    ci_upper = as.numeric(att + 1.96 * se),
    p_value = as.numeric(2 * pnorm(-abs(att / se))),
    n_obs = sum(!is.na(Y_wide)),
    weights = list(theta = theta, omega = omega),
    lambda = list(time = lambda_time, unit = lambda_unit, nn = lambda_nn)
  )
}

# =============================================================================
# RUN ALL ESTIMATORS
# =============================================================================

run_all_estimators <- function(panel, outcome_col = "Y_any", true_att = NULL) {
  results <- list(
    twfe = tryCatch(
      estimate_twfe(panel, outcome_col = outcome_col),
      error = function(e) {
        message("TWFE failed: ", e$message)
        list(method = "TWFE", estimate = NA_real_, se = NA_real_,
             p_value = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
             n_obs = NA_integer_)
      }
    ),
    mc = tryCatch(
      estimate_mc(panel, outcome_col = outcome_col),
      error = function(e) {
        message("MC failed: ", e$message)
        list(method = "Matrix Completion", estimate = NA_real_, se = NA_real_,
             p_value = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
             n_obs = NA_integer_)
      }
    ),
    sdid = tryCatch(
      estimate_sdid(panel, outcome_col = outcome_col),
      error = function(e) {
        message("SDID failed: ", e$message)
        list(method = "Synthetic DiD", estimate = NA_real_, se = NA_real_,
             p_value = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
             n_obs = NA_integer_)
      }
    ),
    trop = tryCatch(
      estimate_trop(panel, outcome_col = outcome_col),
      error = function(e) {
        message("TROP failed: ", e$message)
        list(method = "TROP", estimate = NA_real_, se = NA_real_,
             p_value = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
             n_obs = NA_integer_)
      }
    )
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
