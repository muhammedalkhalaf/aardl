#' Augmented ARDL Cointegration Analysis
#'
#' Estimates an Augmented Autoregressive Distributed Lag (A-ARDL) model and
#' tests for cointegration using the three-test framework of Sam, McNown and
#' Goh (2019). Eight model variants are supported, covering standard ARDL,
#' bootstrap ARDL, Fourier ARDL, and their nonlinear (NARDL) counterparts.
#'
#' @param formula A formula of the form \code{y ~ x1 + x2 + ...}. The
#'   dependent variable is the first variable on the left-hand side.
#' @param data A data frame containing all variables in \code{formula}.
#' @param type Character string specifying the model variant. One of
#'   \code{"aardl"} (default), \code{"baardl"}, \code{"faardl"},
#'   \code{"fbaardl"}, \code{"nardl"}, \code{"fanardl"}, \code{"banardl"},
#'   \code{"fbanardl"}.
#' @param decompose Character vector naming variables to decompose into
#'   positive and negative partial sums (required for NARDL variants).
#' @param max_lag Positive integer. Maximum lag order for selection. Default
#'   \code{4}.
#' @param max_k Numeric. Maximum Fourier frequency for selection in Fourier
#'   variants. Default \code{3}.
#' @param ic Character. Information criterion: \code{"bic"} (default) or
#'   \code{"aic"}.
#' @param reps Positive integer. Number of bootstrap replications for
#'   bootstrap variants. Default \code{999}.
#' @param level Numeric (0, 100). Confidence level. Default \code{95}.
#' @param case Integer 1--5. PSS deterministic case. Default \code{3}
#'   (unrestricted constant, no trend).
#' @param bootstrap Character. Bootstrap method for bootstrap variants:
#'   \code{"bvz"} (Bertelli, Vacca & Zoia 2022, default) or
#'   \code{"mcnown"} (McNown, Sam & Goh 2018).
#'
#' @return An object of class \code{"aardl"} with components:
#'   \describe{
#'     \item{coefficients}{Named vector of EC-form coefficients (ADJ, LR, SR).}
#'     \item{vcov}{Variance-covariance matrix of EC-form coefficients.}
#'     \item{F_pss}{F-statistic: joint significance of all lagged levels.}
#'     \item{t_pss}{t-statistic on the lagged dependent variable level.}
#'     \item{F_ind}{F-statistic: joint significance of lagged independent
#'       variable levels (the Sam 2019 augmentation).}
#'     \item{p_F_pss}{p-value of \code{F_pss}.}
#'     \item{p_t_pss}{Two-sided p-value of \code{t_pss}.}
#'     \item{p_F_ind}{p-value of \code{F_ind}.}
#'     \item{coint_status}{Character: one of \code{"cointegrated"},
#'       \code{"degenerate_1"}, \code{"degenerate_2"},
#'       \code{"no_cointegration"}.}
#'     \item{r2}{R-squared of the EC regression.}
#'     \item{r2_adj}{Adjusted R-squared.}
#'     \item{aic}{AIC of the EC regression.}
#'     \item{bic}{BIC of the EC regression.}
#'     \item{nobs}{Number of observations.}
#'     \item{opt_p}{Selected AR lag for the dependent variable.}
#'     \item{opt_q}{Named integer vector of selected lags for each regressor.}
#'     \item{kstar}{Optimal Fourier frequency (Fourier models only).}
#'     \item{type}{Model type.}
#'     \item{depvar}{Name of the dependent variable.}
#'     \item{indepvars}{Names of the independent variables.}
#'     \item{ic}{Information criterion used.}
#'     \item{bootstrap_method}{Bootstrap method used (bootstrap models only).}
#'   }
#'
#' @references
#' Sam, C. Y., McNown, R. and Goh, S. K. (2019). An augmented autoregressive
#' distributed lag bounds test for cointegration. \emph{Economics Letters},
#' 174, 47--50. \doi{10.1016/j.econlet.2018.12.007}
#'
#' McNown, R., Sam, C. Y. and Goh, S. K. (2018). Bootstrapping the Autoregressive
#' Distributed Lag Test for Cointegration. \emph{Applied Economics}, 50(13),
#' 1509--1521. \doi{10.1080/00036846.2017.1366150}
#'
#' Bertelli, S., Vacca, G. and Zoia, M. G. (2022). Bootstrapping the ARDL
#' Bounds Cointegration Test: A Conditional Approach. \emph{Economics Letters},
#' 216, 110662. \doi{10.1016/j.econlet.2022.110662}
#'
#' Yilanci, V., Bozoklu, S. and Gorus, M. S. (2020). Are OECD Countries
#' Converging in Per Capita Healthcare Expenditures? Evidence from a Fourier
#' LM Unit Root Test. \emph{Emerging Markets Finance and Trade}, 56(3), 642--655.
#' \doi{10.1080/1540496X.2019.1594204}
#'
#' Shin, Y., Yu, B. and Greenwood-Nimmo, M. (2014). Modelling Asymmetric
#' Cointegration and Dynamic Multipliers in a Nonlinear ARDL Framework.
#' In Horrace, W. C. and Sickles, R. C. (Eds.), \emph{Festschrift in Honor of
#' Peter Schmidt}. Springer, New York. \doi{10.1007/978-1-4899-8008-3_9}
#'
#' Pesaran, M. H., Shin, Y. and Smith, R. J. (2001). Bounds testing approaches
#' to the analysis of level relationships. \emph{Journal of Applied Econometrics},
#' 16(3), 289--326. \doi{10.1002/jae.616}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n  <- 80
#' t  <- seq_len(n)
#' x1 <- cumsum(rnorm(n))
#' x2 <- cumsum(rnorm(n))
#' y  <- 0.5 * x1 + 0.3 * x2 + rnorm(n, sd = 0.5)
#' df <- data.frame(y = y, x1 = x1, x2 = x2)
#' result <- aardl(y ~ x1 + x2, data = df, max_lag = 2, ic = "aic")
#' print(result)
#' }
#'
#' @importFrom stats lm coef vcov AIC BIC logLik residuals fitted
#'   reformulate model.matrix model.frame pf pt df.residual
#'   na.omit setNames rnorm
#' @export
aardl <- function(formula, data, type = "aardl",
                  decompose = NULL, max_lag = 4L, max_k = 3,
                  ic = "bic", reps = 999L, level = 95,
                  case = 3L, bootstrap = "bvz") {

  ## --- Input validation ---
  type      <- tolower(as.character(type))
  ic        <- tolower(as.character(ic))
  bootstrap <- tolower(as.character(bootstrap))

  valid_types <- c("aardl", "baardl", "faardl", "fbaardl",
                   "nardl", "fanardl", "banardl", "fbanardl")
  if (!type %in% valid_types) {
    stop(paste("'type' must be one of:", paste(valid_types, collapse = ", ")),
         call. = FALSE)
  }
  if (!ic %in% c("aic", "bic")) {
    stop("'ic' must be \"aic\" or \"bic\".", call. = FALSE)
  }
  if (!bootstrap %in% c("bvz", "mcnown")) {
    stop("'bootstrap' must be \"bvz\" or \"mcnown\".", call. = FALSE)
  }
  if (!case %in% 1L:5L) {
    stop("'case' must be an integer from 1 to 5.", call. = FALSE)
  }

  is_nardl     <- type %in% c("nardl", "fanardl", "banardl", "fbanardl")
  has_bootstrap <- type %in% c("baardl", "fbaardl", "banardl", "fbanardl")
  has_fourier  <- type %in% c("faardl", "fbaardl", "fanardl", "fbanardl")

  if (is_nardl && (is.null(decompose) || length(decompose) == 0L)) {
    stop("NARDL types require 'decompose' to be specified.", call. = FALSE)
  }
  if (!is_nardl && !is.null(decompose)) {
    stop("'decompose' is only for NARDL types.", call. = FALSE)
  }

  ## --- Parse formula ---
  mf       <- model.frame(formula, data = data, na.action = na.omit)
  depvar   <- names(mf)[1L]
  indepvars <- names(mf)[-1L]
  y        <- as.numeric(mf[[depvar]])
  n_full   <- nrow(mf)

  if (n_full < 30L) {
    stop(paste("Too few observations (", n_full, "). Need at least 30."),
         call. = FALSE)
  }

  ## --- NARDL decomposition ---
  indepvars_all <- indepvars  # full list (replaced for NARDL)
  xmat_levels <- as.matrix(mf[, indepvars, drop = FALSE])

  if (is_nardl) {
    pos_neg_list <- .aardl_decompose(mf, decompose)
    xmat_levels  <- pos_neg_list$xmat  # replaces decomposed vars with +/-
    indepvars_all <- pos_neg_list$names
  }
  n_x <- ncol(xmat_levels)

  ## --- Fourier frequency selection ---
  kstar <- 0
  if (has_fourier) {
    kstar <- .aardl_select_fourier(y, xmat_levels, max_k = max_k,
                                   max_lag = max_lag)
  }
  fourier_mat <- if (has_fourier && kstar > 0) {
    .aardl_fourier_terms(kstar = kstar, n = n_full)
  } else {
    NULL
  }

  ## --- Lag selection ---
  lag_result <- .aardl_select_lags(
    y = y, xmat = xmat_levels, max_lag = max_lag,
    ic = ic, fourier_mat = fourier_mat
  )
  opt_p <- lag_result$opt_p
  opt_q <- lag_result$opt_q  # length n_x vector

  ## --- Final regression ---
  reg <- .aardl_build_ecm(
    y = y, xmat = xmat_levels, p = opt_p, q_vec = opt_q,
    fourier_mat = fourier_mat
  )
  fit <- lm(reg$y_dep ~ reg$X_mat - 1)
  n_used  <- length(reg$y_dep)
  df_res  <- fit$df.residual
  res_fit <- residuals(fit)

  ## --- Cointegration test statistics ---
  b_ols <- coef(fit)
  ## F_pss: all lagged level columns
  lv_cols <- reg$level_cols  # indices of lagged level columns in X_mat
  ## t_pss: coefficient on lagged y
  ydv_col <- reg$ydv_col
  ## F_ind: lagged independent level columns
  xlv_cols <- reg$xlv_cols

  ## F_pss (joint test of all lagged levels)
  F_pss_res <- .aardl_ftest(fit, lv_cols)
  F_pss   <- F_pss_res$F
  p_F_pss <- F_pss_res$p

  ## t_pss (t-stat on lagged dependent variable)
  t_pss   <- summary(fit)$coefficients[ydv_col, "t value"]
  p_t_pss <- 2 * pt(-abs(t_pss), df = df_res)

  ## F_ind (joint test of lagged independent levels)
  F_ind_res <- .aardl_ftest(fit, xlv_cols)
  F_ind   <- F_ind_res$F
  p_F_ind <- F_ind_res$p

  ## --- Long-run (EC-form) coefficients via nlcom ---
  adj_coef  <- b_ols[[ydv_col]]
  lr_coefs  <- vapply(xlv_cols, function(j) {
    -b_ols[[j]] / adj_coef
  }, numeric(1L))
  names(lr_coefs) <- paste0("LR.", indepvars_all)

  ## --- Information criteria ---
  ll_val  <- as.numeric(logLik(fit))
  k_fit   <- length(b_ols)
  aic_val <- -2 * ll_val + 2 * k_fit
  bic_val <- -2 * ll_val + k_fit * log(n_used)
  r2      <- summary(fit)$r.squared
  r2_adj  <- summary(fit)$adj.r.squared

  ## --- Cointegration conclusion ---
  sig_fov  <- p_F_pss < (1 - level / 100)
  sig_tdv  <- p_t_pss < (1 - level / 100)
  sig_find <- p_F_ind < (1 - level / 100)

  coint_status <- if (sig_fov && sig_tdv && sig_find) {
    "cointegrated"
  } else if (sig_fov && sig_tdv && !sig_find) {
    "degenerate_2"
  } else if (sig_fov && !sig_tdv && sig_find) {
    "degenerate_1"
  } else {
    "no_cointegration"
  }

  ## --- NARDL asymmetry tests ---
  asym_tests <- NULL
  if (is_nardl) {
    asym_tests <- .aardl_asymmetry_tests(fit, reg, decompose)
  }

  ## --- Bootstrap (if requested) ---
  bs_results <- NULL
  if (has_bootstrap) {
    bs_results <- .aardl_bootstrap(fit, reg, reps = reps,
                                   method = bootstrap, level = level)
    if (!is.null(bs_results)) {
      p_F_pss <- bs_results$p_F_pss
      p_t_pss <- bs_results$p_t_pss
      p_F_ind <- bs_results$p_F_ind
      sig_fov  <- p_F_pss < (1 - level / 100)
      sig_tdv  <- p_t_pss < (1 - level / 100)
      sig_find <- p_F_ind < (1 - level / 100)
      coint_status <- if (sig_fov && sig_tdv && sig_find) {
        "cointegrated"
      } else if (sig_fov && sig_tdv && !sig_find) {
        "degenerate_2"
      } else if (sig_fov && !sig_tdv && sig_find) {
        "degenerate_1"
      } else {
        "no_cointegration"
      }
    }
  }

  ## --- Return ---
  result <- list(
    coefficients    = c(ADJ = adj_coef, lr_coefs),
    vcov            = vcov(fit),
    F_pss           = F_pss,
    t_pss           = t_pss,
    F_ind           = F_ind,
    p_F_pss         = p_F_pss,
    p_t_pss         = p_t_pss,
    p_F_ind         = p_F_ind,
    coint_status    = coint_status,
    r2              = r2,
    r2_adj          = r2_adj,
    aic             = aic_val,
    bic             = bic_val,
    nobs            = n_used,
    opt_p           = opt_p,
    opt_q           = setNames(opt_q, indepvars_all),
    kstar           = kstar,
    type            = type,
    depvar          = depvar,
    indepvars       = indepvars,
    indepvars_all   = indepvars_all,
    ic              = ic,
    bootstrap_method = if (has_bootstrap) bootstrap else NULL,
    bootstrap_reps   = if (has_bootstrap) reps else NULL,
    bs_results       = bs_results,
    asym_tests       = asym_tests,
    lm_fit           = fit
  )
  class(result) <- "aardl"
  result
}


## ===========================================================================
## INTERNAL HELPERS
## ===========================================================================

#' NARDL partial sum decomposition
#' @keywords internal
.aardl_decompose <- function(mf, decompose) {
  n <- nrow(mf)
  new_cols <- list()
  new_names <- character(0L)

  indepvars <- names(mf)[-1L]

  for (xv in indepvars) {
    if (xv %in% decompose) {
      dx   <- c(NA_real_, diff(as.numeric(mf[[xv]])))
      xpos <- cumsum(ifelse(is.na(dx) | dx < 0, 0, dx))
      xneg <- cumsum(ifelse(is.na(dx) | dx > 0, 0, dx))
      pname <- paste0(xv, "_pos")
      nname <- paste0(xv, "_neg")
      new_cols[[pname]] <- xpos
      new_cols[[nname]] <- xneg
      new_names <- c(new_names, pname, nname)
    } else {
      new_cols[[xv]] <- as.numeric(mf[[xv]])
      new_names <- c(new_names, xv)
    }
  }
  list(xmat = do.call(cbind, new_cols), names = new_names)
}


#' Fourier terms matrix
#' @keywords internal
.aardl_fourier_terms <- function(kstar, n) {
  t_idx <- seq_len(n)
  sin_t <- sin(2 * pi * kstar * t_idx / n)
  cos_t <- cos(2 * pi * kstar * t_idx / n)
  cbind(sin_ft = sin_t, cos_ft = cos_t)
}


#' Select Fourier frequency by minimum SSR
#' @keywords internal
.aardl_select_fourier <- function(y, xmat, max_k, max_lag) {
  n    <- length(y)
  k_grid <- seq(0.1, max_k, by = 0.1)
  best_ssr <- Inf
  best_k   <- k_grid[1L]

  for (kv in k_grid) {
    fm <- .aardl_fourier_terms(kstar = kv, n = n)
    reg <- .aardl_build_ecm(y, xmat, p = 1L,
                             q_vec = rep(1L, ncol(xmat)),
                             fourier_mat = fm)
    fit <- tryCatch(lm(reg$y_dep ~ reg$X_mat - 1),
                    error = function(e) NULL)
    if (is.null(fit)) next
    ssr <- sum(residuals(fit)^2)
    if (ssr < best_ssr) {
      best_ssr <- ssr
      best_k   <- kv
    }
  }
  best_k
}


#' Select lag orders by IC
#' @keywords internal
.aardl_select_lags <- function(y, xmat, max_lag, ic, fourier_mat) {
  n_x    <- ncol(xmat)
  best_ic <- Inf
  opt_p   <- 1L
  opt_q   <- rep(1L, n_x)

  for (p in 1L:max_lag) {
    q_combos <- expand.grid(lapply(seq_len(n_x), function(i) 0L:max_lag))
    for (ci in seq_len(nrow(q_combos))) {
      q_vec <- as.integer(q_combos[ci, ])
      reg   <- .aardl_build_ecm(y, xmat, p = p, q_vec = q_vec,
                                 fourier_mat = fourier_mat)
      if (nrow(reg$X_mat) < ncol(reg$X_mat) + 10L) next
      fit <- tryCatch(lm(reg$y_dep ~ reg$X_mat - 1),
                      error = function(e) NULL)
      if (is.null(fit)) next
      n_u  <- length(reg$y_dep)
      k_u  <- length(coef(fit))
      ll_u <- as.numeric(logLik(fit))
      ic_v <- if (ic == "aic") {
        -2 * ll_u + 2 * k_u
      } else {
        -2 * ll_u + k_u * log(n_u)
      }
      if (ic_v < best_ic) {
        best_ic <- ic_v
        opt_p   <- p
        opt_q   <- q_vec
      }
    }
  }
  list(opt_p = opt_p, opt_q = opt_q)
}


#' Build ECM regression matrices
#'
#' Returns lagged-level and differenced matrices, plus column index vectors.
#' @keywords internal
.aardl_build_ecm <- function(y, xmat, p, q_vec, fourier_mat = NULL) {
  n    <- length(y)
  n_x  <- ncol(xmat)
  ## Maximum lag needed
  max_lag_used <- max(p, max(q_vec, 0L)) + 1L
  if (n <= max_lag_used + 5L) {
    return(list(y_dep = numeric(0), X_mat = matrix(nrow = 0, ncol = 0),
                level_cols = integer(0), ydv_col = 1L, xlv_cols = integer(0)))
  }

  dy <- diff(y)  # length n-1

  ## Number of rows after differencing and lags
  n_trim <- n - max_lag_used
  idx    <- (max_lag_used + 1L):(n - 1L)  # indices in dy
  if (length(idx) < 5L) {
    return(list(y_dep = numeric(0), X_mat = matrix(nrow = 0, ncol = 0),
                level_cols = integer(0), ydv_col = 1L, xlv_cols = integer(0)))
  }
  n_used <- length(idx)

  y_dep <- dy[idx]

  ## Lagged level of y: y_{t-1}
  y_lag1 <- y[idx]   # y at time t-1 (dy[idx] = y[idx+1] - y[idx])
  ## Constant
  intercept <- rep(1, n_used)

  ## Lagged levels of x
  x_lag1_mat <- xmat[idx, , drop = FALSE]

  ## Lagged differences of y (SR): L1Dy, ..., Lp.Dy
  ly_diff_mat <- matrix(NA_real_, nrow = n_used, ncol = p)
  for (j in seq_len(p)) {
    ly_diff_mat[, j] <- dy[idx - j]
  }
  colnames(ly_diff_mat) <- paste0("Ldy_", seq_len(p))

  ## Contemp + lagged diffs of x
  x_diff_list <- list()
  for (i in seq_len(n_x)) {
    dx_i <- diff(xmat[, i])
    qi   <- q_vec[i]
    for (j in 0L:qi) {
      nm <- if (j == 0L) paste0("Dx_", i) else paste0("Ldy_x", i, "_", j)
      x_diff_list[[nm]] <- dx_i[idx - j]
    }
  }
  x_diff_mat <- do.call(cbind, x_diff_list)

  ## Track column indices
  ## Column layout: intercept | y_lag1 | x_lag1_mat (n_x cols) |
  ##                ly_diff_mat (p cols) | x_diff_mat | fourier
  col_intercept <- 1L
  col_ydv       <- 2L
  col_xlv_start <- 3L
  col_xlv_end   <- 2L + n_x

  all_mats <- list(intercept = matrix(intercept, ncol = 1L,
                                      dimnames = list(NULL, "intercept")),
                   y_lag1    = matrix(y_lag1, ncol = 1L,
                                      dimnames = list(NULL, "y_lag1")),
                   x_lag1    = x_lag1_mat,
                   ly_diff   = ly_diff_mat,
                   x_diff    = x_diff_mat)

  if (!is.null(fourier_mat)) {
    fm_sub <- fourier_mat[idx, , drop = FALSE]
    all_mats$fourier <- fm_sub
  }

  X_mat <- do.call(cbind, all_mats)

  level_cols <- c(col_ydv, seq(col_xlv_start, col_xlv_end))
  xlv_cols   <- seq(col_xlv_start, col_xlv_end)

  list(y_dep     = y_dep,
       X_mat     = X_mat,
       level_cols = level_cols,
       ydv_col   = col_ydv,
       xlv_cols  = xlv_cols)
}


#' F-test for a subset of coefficients
#' @keywords internal
.aardl_ftest <- function(fit, col_indices) {
  if (length(col_indices) == 0L) {
    return(list(F = NA_real_, p = NA_real_))
  }
  ## Restricted model: set selected coefficients to 0
  b   <- coef(fit)
  V   <- vcov(fit)
  idx <- col_indices
  R   <- matrix(0, nrow = length(idx), ncol = length(b))
  for (i in seq_along(idx)) R[i, idx[i]] <- 1
  r   <- rep(0, length(idx))
  Rb  <- R %*% b - r
  RVR <- R %*% V %*% t(R)
  tryCatch({
    F_stat <- as.numeric(t(Rb) %*% solve(RVR) %*% Rb) / length(idx)
    df_r   <- fit$df.residual
    p_val  <- pf(F_stat, df1 = length(idx), df2 = df_r,
                 lower.tail = FALSE)
    list(F = F_stat, p = p_val)
  }, error = function(e) {
    list(F = NA_real_, p = NA_real_)
  })
}


#' NARDL asymmetry tests (Wald)
#' @keywords internal
.aardl_asymmetry_tests <- function(fit, reg, decompose) {
  results <- list()
  b <- coef(fit)
  V <- vcov(fit)
  nm <- colnames(reg$X_mat)

  for (xv in decompose) {
    pname_lr <- paste0(xv, "_pos")
    nname_lr <- paste0(xv, "_neg")
    p_col <- which(nm == pname_lr)
    n_col <- which(nm == nname_lr)
    if (length(p_col) == 1L && length(n_col) == 1L) {
      diff_b  <- b[p_col] - b[n_col]
      diff_var <- V[p_col, p_col] + V[n_col, n_col] - 2 * V[p_col, n_col]
      W <- diff_b^2 / diff_var
      p_w <- pf(W, df1 = 1L, df2 = fit$df.residual, lower.tail = FALSE)
      results[[paste0(xv, "_LR")]] <- list(W = W, p = p_w)
    }
    dp_nm <- paste0("Dx_", which(colnames(reg$xmat_levels) ==
                                    pname_lr))
    dn_nm <- paste0("Dx_", which(colnames(reg$xmat_levels) ==
                                    nname_lr))
    dp_col <- which(nm == dp_nm)
    dn_col <- which(nm == dn_nm)
    if (length(dp_col) == 1L && length(dn_col) == 1L) {
      diff_b  <- b[dp_col] - b[dn_col]
      diff_var <- V[dp_col, dp_col] + V[dn_col, dn_col] -
                  2 * V[dp_col, dn_col]
      W <- diff_b^2 / diff_var
      p_w <- pf(W, df1 = 1L, df2 = fit$df.residual, lower.tail = FALSE)
      results[[paste0(xv, "_SR")]] <- list(W = W, p = p_w)
    }
  }
  results
}


#' Parametric bootstrap for ARDL bounds tests
#' @keywords internal
.aardl_bootstrap <- function(fit, reg, reps, method, level) {
  y_dep <- reg$y_dep
  X_mat <- reg$X_mat
  n_used <- length(y_dep)
  b_ols  <- coef(fit)
  res    <- residuals(fit)

  F_pss_vec <- numeric(reps)
  t_pss_vec <- numeric(reps)
  F_ind_vec <- numeric(reps)

  for (r in seq_len(reps)) {
    if (method == "bvz") {
      ## Conditional bootstrap: resample residuals
      res_b <- sample(res, size = n_used, replace = TRUE)
    } else {
      ## McNown: unconditional normal bootstrap
      res_b <- rnorm(n_used, sd = sqrt(mean(res^2)))
    }
    y_b <- as.numeric(X_mat %*% b_ols) + res_b
    fit_b <- tryCatch(lm(y_b ~ X_mat - 1), error = function(e) NULL)
    if (is.null(fit_b)) next
    F_pss_vec[r] <- .aardl_ftest(fit_b, reg$level_cols)$F
    t_pss_vec[r] <- tryCatch(
      summary(fit_b)$coefficients[reg$ydv_col, "t value"],
      error = function(e) NA_real_
    )
    F_ind_vec[r] <- .aardl_ftest(fit_b, reg$xlv_cols)$F
  }

  obs_F_pss <- .aardl_ftest(fit, reg$level_cols)$F
  obs_t_pss <- summary(fit)$coefficients[reg$ydv_col, "t value"]
  obs_F_ind <- .aardl_ftest(fit, reg$xlv_cols)$F

  p_F_pss <- mean(F_pss_vec >= obs_F_pss, na.rm = TRUE)
  p_t_pss <- mean(abs(t_pss_vec) >= abs(obs_t_pss), na.rm = TRUE)
  p_F_ind <- mean(F_ind_vec >= obs_F_ind, na.rm = TRUE)

  list(p_F_pss = p_F_pss, p_t_pss = p_t_pss, p_F_ind = p_F_ind,
       reps = reps)
}
