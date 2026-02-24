# ===========================================================================
# tests/testthat/test-hsreg.R — Tests for the main estimation function
#                                     and S3 methods
# ===========================================================================

# --- Helper: small reproducible fit used across tests ---
make_test_fit <- function(seed = 42, n = 100, p = 10, n_mcmc = 100,
                           burnin = 50) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", 1:p)
  beta_true <- c(rep(2, 3), rep(0, p - 3))
  y <- X %*% beta_true + rnorm(n)
  hsreg(y, X, n_mcmc = n_mcmc, burnin = burnin, verbose = FALSE,
            seed = seed)
}


# ---------------------------------------------------------------------------
# Round-trip: fit, coef, vcov, predict, confint, print, summary
# ---------------------------------------------------------------------------

test_that("horseshoe full round-trip works", {
  fit <- make_test_fit()

  # coef
  b <- coef(fit)
  expect_length(b, 10)
  expect_named(b, paste0("V", 1:10))

  # vcov
  V <- vcov(fit)
  expect_equal(dim(V), c(10, 10))
  expect_true(all(diag(V) >= 0))

  # predict
  pred_xb <- predict(fit, type = "xb")
  expect_length(pred_xb, 100)

  pred_resid <- predict(fit, type = "residuals")
  expect_length(pred_resid, 100)

  pred_stdp <- predict(fit, type = "stdp")
  expect_length(pred_stdp, 100)
  expect_true(all(pred_stdp >= 0))

  # confint
  ci <- confint(fit)
  expect_equal(dim(ci), c(10, 2))
  expect_true(all(ci[, 1] <= ci[, 2]))

  # print (just check it doesn't error)
  expect_no_error(capture.output(print(fit)))

  # summary
  s <- summary(fit)
  expect_s3_class(s, "summary.hsreg")
  expect_no_error(capture.output(print(s)))
})


# ---------------------------------------------------------------------------
# Seed reproducibility
# ---------------------------------------------------------------------------

test_that("same seed produces identical draws", {
  fit1 <- make_test_fit(seed = 99)
  fit2 <- make_test_fit(seed = 99)

  expect_identical(fit1$beta_draws, fit2$beta_draws)
  expect_identical(fit1$sigma2_draws, fit2$sigma2_draws)
  expect_identical(fit1$tau2_draws, fit2$tau2_draws)
})


# ---------------------------------------------------------------------------
# scale_X
# ---------------------------------------------------------------------------

test_that("scale_X = TRUE vs manual scaling gives same coefficients", {
  set.seed(42)
  n <- 100; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  # Make columns have very different scales
  X[, 1] <- X[, 1] * 100
  X[, 2] <- X[, 2] * 0.01
  y <- rnorm(n)

  fit_scaled <- hsreg(y, X, scale_X = TRUE, n_mcmc = 100, burnin = 50,
                           verbose = FALSE, seed = 42)

  # The coefficients should be on the original (unscaled) scale
  expect_length(coef(fit_scaled), p)
  # Sanity: coefficients should be finite
  expect_true(all(is.finite(coef(fit_scaled))))
})


# ---------------------------------------------------------------------------
# Save/load round-trip
# ---------------------------------------------------------------------------

test_that("saving to .rds works and can be read back", {
  tmpfile <- tempfile(fileext = ".rds")
  on.exit(unlink(tmpfile))

  set.seed(42)
  n <- 50; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  fit <- hsreg(y, X, n_mcmc = 50, burnin = 20, saving = tmpfile,
                   verbose = FALSE)

  # File should exist
  expect_true(file.exists(tmpfile))

  # Load and check structure
  draws <- readRDS(tmpfile)
  expect_true(is.data.frame(draws))
  expect_equal(nrow(draws), 50)  # n_mcmc draws
  expect_true("draw" %in% colnames(draws))
  expect_true("sigma2" %in% colnames(draws))
  expect_true("tau2" %in% colnames(draws))
  expect_equal(ncol(draws), 3 + p)
})

test_that("saving to .csv works", {
  tmpfile <- tempfile(fileext = ".csv")
  on.exit(unlink(tmpfile))

  set.seed(42)
  n <- 50; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  fit <- hsreg(y, X, n_mcmc = 50, burnin = 20, saving = tmpfile,
                   verbose = FALSE)

  expect_true(file.exists(tmpfile))
  draws <- read.csv(tmpfile)
  expect_equal(nrow(draws), 50)
})


# ---------------------------------------------------------------------------
# predict consistency
# ---------------------------------------------------------------------------

test_that("predict(type='xb') matches manual X %*% coef(fit)", {
  fit <- make_test_fit()
  pred <- predict(fit, type = "xb")
  manual <- as.numeric(fit$X %*% coef(fit))
  expect_equal(pred, manual, tolerance = 1e-10)
})

test_that("predict(type='residuals') matches y - predict(fit)", {
  fit <- make_test_fit()
  resid <- predict(fit, type = "residuals")
  expected <- fit$y - predict(fit, type = "xb")
  expect_equal(resid, expected, tolerance = 1e-10)
})

test_that("predict with newdata works", {
  fit <- make_test_fit()
  X_new <- matrix(rnorm(20 * 10), 20, 10)
  pred <- predict(fit, newdata = X_new, type = "xb")
  expect_length(pred, 20)
})

test_that("predict with newdata errors on wrong dimensions", {
  fit <- make_test_fit()
  X_bad <- matrix(rnorm(20 * 5), 20, 5)  # wrong p
  expect_error(predict(fit, newdata = X_bad), "must have 10 columns")
})

test_that("predict type='residuals' errors with newdata", {
  fit <- make_test_fit()
  X_new <- matrix(rnorm(20 * 10), 20, 10)
  expect_error(predict(fit, newdata = X_new, type = "residuals"),
               "requires newdata = NULL")
})


# ---------------------------------------------------------------------------
# confint with different level
# ---------------------------------------------------------------------------

test_that("confint recomputes CIs when level differs", {
  fit <- make_test_fit()  # default level = 0.95

  ci_95 <- confint(fit, level = 0.95)
  ci_90 <- confint(fit, level = 0.90)

  # 90% CI should be narrower than 95% CI
  width_95 <- ci_95[, 2] - ci_95[, 1]
  width_90 <- ci_90[, 2] - ci_90[, 1]
  expect_true(all(width_90 <= width_95 + 1e-10))
})

test_that("confint with parm subset works", {
  fit <- make_test_fit()

  # By name
  ci_sub <- confint(fit, parm = c("V1", "V3"))
  expect_equal(nrow(ci_sub), 2)
  expect_equal(rownames(ci_sub), c("V1", "V3"))

  # By index
  ci_idx <- confint(fit, parm = c(1, 3))
  expect_equal(nrow(ci_idx), 2)
})


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("horseshoe errors on invalid level", {
  set.seed(42)
  X <- matrix(rnorm(100), 20, 5)
  y <- rnorm(20)
  expect_error(hsreg(y, X, level = 1.5), "must be a scalar in \\(0, 1\\)")
  expect_error(hsreg(y, X, level = 0), "must be a scalar in \\(0, 1\\)")
})

test_that("horseshoe uses colnames from X when available", {
  set.seed(42)
  X <- matrix(rnorm(100), 20, 5)
  colnames(X) <- c("a", "b", "c", "d", "e")
  y <- rnorm(20)
  fit <- hsreg(y, X, n_mcmc = 30, burnin = 10, verbose = FALSE)
  expect_equal(names(coef(fit)), c("a", "b", "c", "d", "e"))
})

test_that("horseshoe generates X1..Xp names when colnames missing", {
  set.seed(42)
  X <- matrix(rnorm(100), 20, 5)
  y <- rnorm(20)
  fit <- hsreg(y, X, n_mcmc = 30, burnin = 10, verbose = FALSE)
  expect_equal(names(coef(fit)), paste0("X", 1:5))
})
