#' ============================================================================
#' Estimators: TWFE, Matrix Completion (fect), Synthetic DiD
#'
#' Target estimand:
#'   By default (treat_type = "post"):
#'     D_it = near_libya_i × post_treat_t
#'   For a short-lived effect window (treat_type = "window"):
#'     D_it = near_libya_i × effect_window_t
#'   ATT on a chosen outcome Y (default: Y_any - binary indicator)
#'
#' Notes:
#'   1) enforce (unit_id, time_id) ordering
#'   2) for wide matrices (SDID) enforce numeric time ordering of columns
#'   3) for fect we allow small lambda so it doesn't collapse to TWFE
#' ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

FECT_AVAILABLE <- requireNamespace("fect", quietly = TRUE)
if (FECT_AVAILABLE) suppressPackageStartupMessages(library(fect))

SYNTHDID_AVAILABLE <- requireNamespace("synthdid", quietly = TRUE)
if (SYNTHDID_AVAILABLE) suppressPackageStartupMessages(library(synthdid))

`%||%` <- function(x, y) if (is.null(x)) y else x

# =============================================================================
# HELPERS (wide matrices with correct time ordering)
# =============================================================================
.as_time_ordered_matrix <- function(data, value_col = "Y") {
  stopifnot(all(c("unit_id", "time_id", value_col) %in% names(data)))

  wide <- data %>%
    select(unit_id, time_id, .data[[value_col]]) %>%
    tidyr::pivot_wider(names_from = time_id, values_from = .data[[value_col]]) %>%
    arrange(unit_id)

  mat <- wide %>%
    select(-unit_id) %>%
    as.matrix()

  cn <- colnames(mat)
  ti <- suppressWarnings(as.integer(cn))
  if (any(is.na(ti))) stop("Non-numeric time_id column names after pivot_wider().")

  mat[, order(ti), drop = FALSE]
}

# =============================================================================
# Prepare data for estimation
# =============================================================================
prepare_estimation_data <- function(
  panel,
  outcome_col = "Y_any",
  treat_type = c("post", "window"),
  trim_after_window = (match.arg(treat_type) == "window")
) {
  treat_type <- match.arg(treat_type)
  stopifnot(outcome_col %in% names(panel))

  if (treat_type == "window" && isTRUE(trim_after_window)) {
    if (!("effect_window" %in% names(panel))) {
      stop("treat_type='window' requires panel$effect_window (from generate_time_panel())")
    }
    panel <- panel %>% filter(post_treat == 0 | effect_window == 1)
  }

  post_est <- if (treat_type == "post") panel$post_treat else panel$effect_window

  panel %>%
    arrange(unit_id, time_id) %>%
    mutate(
      post_est = as.integer(post_est),
      D = as.integer(near_libya * post_est),
      Y = as.numeric(.data[[outcome_col]])
    )
}

# =============================================================================
# TWFE
# =============================================================================
estimate_twfe <- function(
  panel,
  outcome_col = "Y_any",
  treat_type = c("post", "window"),
  trim_after_window = NULL,
  cluster = "unit_id"
) {
  treat_type <- match.arg(treat_type)
  if (is.null(trim_after_window)) trim_after_window <- (treat_type == "window")

  data <- prepare_estimation_data(
    panel = panel,
    outcome_col = outcome_col,
    treat_type = treat_type,
    trim_after_window = trim_after_window
  )

  fit <- feols(Y ~ D | unit_id + time_id, data = data, cluster = cluster)

  ct <- as.data.frame(fixest::coeftable(fit))
  est <- as.numeric(ct["D", "Estimate"])
  se <- as.numeric(ct["D", "Std. Error"])
  pval <- ct["D", "Pr(>|t|)"] %||% ct["D", "Pr(>|z|)"]

  list(
    method = "TWFE",
    estimand = outcome_col,
    estimate = est,
    se = se,
    ci_lower = est - 1.96 * se,
    ci_upper = est + 1.96 * se,
    p_value = as.numeric(pval),
    n_obs = fit$nobs,
    fit = fit,
    se_fallback = FALSE
  )
}

# =============================================================================
# MATRIX COMPLETION (fect)
# =============================================================================
.build_lambda_grid <- function(lambda_base = 0.005) {
  base <- as.numeric(lambda_base)[1]
  # Expanded grid: added larger values (1, 2) for strong latent trends

  mult <- c(400, 200, 64, 32, 16, 8, 4, 2, 1, 0.5, 0.25, 0.125, 0.0625, 0.03125)
  sort(unique(pmax(base * mult, 1e-5)))
}

estimate_mc <- function(
  panel,
  outcome_col = "Y_any",
  treat_type = c("post", "window"),
  trim_after_window = NULL,
  lambda = 0.005,
  use_cv = TRUE,
  n_boot = 100,
  k_cv = 5,
  cv_nobs = 1,
  parallel = FALSE,
  seed = NULL,
  quiet = TRUE
) {
  if (!FECT_AVAILABLE) stop("fect package not installed. Run install.packages('fect')")

  treat_type <- match.arg(treat_type)
  if (is.null(trim_after_window)) trim_after_window <- (treat_type == "window")

  data <- prepare_estimation_data(
    panel = panel,
    outcome_col = outcome_col,
    treat_type = treat_type,
    trim_after_window = trim_after_window
  )

  df <- data %>%
    select(unit_id, time_id, Y, D) %>%
    mutate(
      unit_id = as.integer(as.factor(unit_id)),
      time_id = as.integer(as.factor(time_id)),
      D = as.integer(D)
    )

  lambda_arg <- if (isTRUE(use_cv)) {
    if (length(lambda) > 1) as.numeric(lambda) else .build_lambda_grid(lambda)
  } else {
    as.numeric(lambda)[1]
  }

  run_fect <- function() {
    if (!isTRUE(quiet)) {
      return(
        fect::fect(
          Y ~ D,
          data = df,
          index = c("unit_id", "time_id"),
          force = "two-way",
          method = "mc",
          lambda = lambda_arg,
          nlambda = if (length(lambda_arg) > 1) length(lambda_arg) else 8,
          CV = isTRUE(use_cv),
          k = k_cv,
          cv.nobs = cv_nobs,
          cv.donut = 0,
          se = TRUE,
          vartype = "bootstrap",
          nboots = n_boot,
          proportion = 0,
          parallel = parallel,
          seed = seed
        )
      )
    }

    tf_out <- tempfile()
    tf_msg <- tempfile()
    con_out <- file(tf_out, open = "wt")
    con_msg <- file(tf_msg, open = "wt")
    sink(con_out, type = "output")
    sink(con_msg, type = "message")
    on.exit(
      {
        try(sink(type = "message"), silent = TRUE)
        try(sink(type = "output"), silent = TRUE)
        try(close(con_msg), silent = TRUE)
        try(close(con_out), silent = TRUE)
        unlink(c(tf_out, tf_msg))
      },
      add = TRUE
    )

    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data = df,
        index = c("unit_id", "time_id"),
        force = "two-way",
        method = "mc",
        lambda = lambda_arg,
        nlambda = if (length(lambda_arg) > 1) length(lambda_arg) else 8,
        CV = isTRUE(use_cv),
        k = k_cv,
        cv.nobs = cv_nobs,
        cv.donut = 0,
        se = TRUE,
        vartype = "bootstrap",
        nboots = n_boot,
        proportion = 0,
        parallel = parallel,
        seed = seed
      )
    ))
  }

  fit <- run_fect()
  att <- as.numeric(fit$att.avg)

  boot <- if (!is.null(fit$att.avg.boot)) as.numeric(fit$att.avg.boot) else NULL
  se <- if (!is.null(boot)) sd(boot, na.rm = TRUE) else NA_real_

  # Track SE fallback

  se_fallback <- FALSE
  if (is.na(se) || se == 0) {
    se <- abs(att) * 0.1 + 0.001
    se_fallback <- TRUE
  }

  if (!is.null(boot) && !se_fallback) {
    qs <- stats::quantile(boot, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    ci_lower <- as.numeric(qs[1])
    ci_upper <- as.numeric(qs[2])
  } else {
    ci_lower <- att - 1.96 * se
    ci_upper <- att + 1.96 * se
  }

  list(
    method = "Matrix Completion",
    estimand = outcome_col,
    estimate = att,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = as.numeric(2 * pnorm(-abs(att / se))),
    n_obs = nrow(df),
    fit = fit,
    se_fallback = se_fallback,
    diagnostics = list(
      lambda_used = lambda_arg,
      lambda_cv = fit$lambda.cv %||% NA,
      MSPE = fit$MSPE %||% NA
    )
  )
}

# =============================================================================
# SYNTHETIC DID
# =============================================================================
estimate_sdid <- function(
  panel,
  outcome_col = "Y_any",
  treat_type = c("post", "window"),
  trim_after_window = NULL,
  se_method = "placebo"
) {
  if (!SYNTHDID_AVAILABLE) stop("synthdid package not installed. Run install.packages('synthdid')")

  treat_type <- match.arg(treat_type)
  if (is.null(trim_after_window)) trim_after_window <- (treat_type == "window")

  data <- prepare_estimation_data(
    panel = panel,
    outcome_col = outcome_col,
    treat_type = treat_type,
    trim_after_window = trim_after_window
  )

  time_info <- data %>%
    distinct(time_id, post_est) %>%
    arrange(time_id)
  T0 <- sum(time_info$post_est == 0)

  Y_mat <- .as_time_ordered_matrix(data, value_col = "Y")

  unit_treated <- data %>%
    group_by(unit_id) %>%
    summarise(treated = max(near_libya), .groups = "drop") %>%
    arrange(unit_id) %>%
    pull(treated)

  N <- nrow(Y_mat)
  N1 <- sum(unit_treated == 1)
  N0 <- N - N1

  ctrl_idx <- which(unit_treated == 0)
  treat_idx <- which(unit_treated == 1)
  Y_reordered <- Y_mat[c(ctrl_idx, treat_idx), , drop = FALSE]

  fit <- tryCatch(synthdid::synthdid_estimate(Y_reordered, N0, T0), error = function(e) NULL)
  if (is.null(fit)) {
    return(list(
      method = "SynthDiD", estimand = outcome_col,
      estimate = NA_real_, se = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      p_value = NA_real_, n_obs = sum(!is.na(Y_mat)),
      se_fallback = FALSE
    ))
  }

  att <- as.numeric(fit)
  se <- tryCatch(as.numeric(synthdid::synthdid_se(fit, method = se_method)), error = function(e) NA_real_)

  # Track SE fallback
  se_fallback <- FALSE
  if (is.na(se) || se == 0) {
    se <- abs(att) * 0.1 + 0.001
    se_fallback <- TRUE
  }

  list(
    method = "SynthDiD",
    estimand = outcome_col,
    estimate = att,
    se = se,
    ci_lower = att - 1.96 * se,
    ci_upper = att + 1.96 * se,
    p_value = as.numeric(2 * pnorm(-abs(att / se))),
    n_obs = sum(!is.na(Y_mat)),
    fit = fit,
    se_fallback = se_fallback
  )
}

# =============================================================================
# RUN ALL ESTIMATORS (clean set)
# =============================================================================
run_all_estimators <- function(
  panel,
  outcome_col = "Y_any",
  true_att = NULL,
  treat_type = c("post", "window"),
  trim_after_window = NULL,
  estimators = c("twfe", "mc", "sdid")
) {
  treat_type <- match.arg(treat_type)
  if (is.null(trim_after_window)) trim_after_window <- (treat_type == "window")

  results <- list()

  if ("twfe" %in% estimators) {
    results$twfe <- tryCatch(
      estimate_twfe(panel, outcome_col = outcome_col, treat_type = treat_type, trim_after_window = trim_after_window),
      error = function(e) list(method = "TWFE", estimate = NA_real_, se = NA_real_, p_value = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, n_obs = NA_integer_, se_fallback = FALSE)
    )
  }

  if ("mc" %in% estimators) {
    results$mc <- tryCatch(
      estimate_mc(panel, outcome_col = outcome_col, treat_type = treat_type, trim_after_window = trim_after_window),
      error = function(e) list(method = "Matrix Completion", estimate = NA_real_, se = NA_real_, p_value = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, n_obs = NA_integer_, se_fallback = FALSE)
    )
  }

  if ("sdid" %in% estimators) {
    if (!SYNTHDID_AVAILABLE) {
      message("Skipping SynthDiD: package not installed")
    } else {
      results$sdid <- tryCatch(
        estimate_sdid(panel, outcome_col = outcome_col, treat_type = treat_type, trim_after_window = trim_after_window),
        error = function(e) list(method = "SynthDiD", estimate = NA_real_, se = NA_real_, p_value = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, n_obs = NA_integer_, se_fallback = FALSE)
      )
    }
  }

  results_df <- purrr::map_dfr(results, function(r) {
    tibble(
      method = r$method,
      estimand = outcome_col,
      estimate = r$estimate,
      se = r$se,
      p_value = r$p_value,
      ci_lower = r$ci_lower %||% (r$estimate - 1.96 * r$se),
      ci_upper = r$ci_upper %||% (r$estimate + 1.96 * r$se),
      n_obs = r$n_obs %||% NA_integer_,
      se_fallback = r$se_fallback %||% FALSE
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
