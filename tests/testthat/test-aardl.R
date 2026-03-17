test_that("aardl returns correct class", {
  set.seed(1)
  n  <- 60
  x1 <- cumsum(rnorm(n))
  y  <- 0.5 * x1 + rnorm(n, sd = 0.3)
  df <- data.frame(y = y, x1 = x1)
  res <- aardl(y ~ x1, data = df, max_lag = 2, ic = "aic")
  expect_s3_class(res, "aardl")
})

test_that("aardl has required output names", {
  set.seed(2)
  n  <- 60
  x1 <- cumsum(rnorm(n))
  y  <- 0.5 * x1 + rnorm(n, sd = 0.3)
  df <- data.frame(y = y, x1 = x1)
  res <- aardl(y ~ x1, data = df, max_lag = 2)
  expect_true(all(c("F_pss", "t_pss", "F_ind",
                    "p_F_pss", "p_t_pss", "p_F_ind",
                    "coint_status", "r2", "aic", "bic",
                    "nobs", "opt_p", "opt_q") %in% names(res)))
})

test_that("aardl validates type", {
  set.seed(3)
  n  <- 50
  df <- data.frame(y = cumsum(rnorm(n)), x1 = cumsum(rnorm(n)))
  expect_error(aardl(y ~ x1, data = df, type = "badtype"),
               "'type' must be one of")
})

test_that("aardl validates NARDL requires decompose", {
  set.seed(3)
  n  <- 50
  df <- data.frame(y = cumsum(rnorm(n)), x1 = cumsum(rnorm(n)))
  expect_error(aardl(y ~ x1, data = df, type = "nardl"),
               "require 'decompose'")
})

test_that("aardl coint_status is a valid category", {
  set.seed(5)
  n  <- 70
  x1 <- cumsum(rnorm(n))
  y  <- 0.8 * x1 + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x1 = x1)
  res <- aardl(y ~ x1, data = df, max_lag = 2)
  expect_true(res$coint_status %in%
    c("cointegrated", "no_cointegration", "degenerate_1", "degenerate_2"))
})

test_that("print.aardl uses message not cat", {
  set.seed(1)
  n  <- 60
  df <- data.frame(y = cumsum(rnorm(n)), x1 = cumsum(rnorm(n)))
  res <- aardl(y ~ x1, data = df, max_lag = 2)
  expect_message(print(res))
  out <- capture.output(print(res))
  expect_equal(length(out), 0L)
})
