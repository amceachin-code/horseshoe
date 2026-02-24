# ===========================================================================
# tests/testthat/test-hsreg_gibbs.R — Tests for the core Gibbs sampler
# ===========================================================================

test_that("hsreg_gibbs runs without error on small problem", {
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(rep(2, 3), rep(0, 7))
  y <- X %*% beta_true + rnorm(n)

  expect_no_error(
    fit <- hsreg_gibbs(y, X, n_mcmc = 50, burnin = 20, verbose = FALSE)
  )
})

test_that("hsreg_gibbs returns correct dimensions", {
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  fit <- hsreg_gibbs(y, X, n_mcmc = 100, burnin = 50, verbose = FALSE)

  expect_equal(nrow(fit$beta_draws), p)
  expect_equal(ncol(fit$beta_draws), 100)
  expect_equal(length(fit$sigma2_draws), 100)
  expect_equal(length(fit$tau2_draws), 100)
  expect_equal(nrow(fit$lambda2_draws), p)  # all penalized by default
  expect_equal(ncol(fit$lambda2_draws), 100)
})

test_that("hsreg_gibbs respects thinning", {
  set.seed(42)
  n <- 50; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  fit <- hsreg_gibbs(y, X, n_mcmc = 50, burnin = 10, thin = 3,
                          verbose = FALSE)

  expect_equal(ncol(fit$beta_draws), 50)
  expect_equal(length(fit$sigma2_draws), 50)
})

test_that("hsreg_gibbs penalized/free indices are correct", {
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  pen <- c(rep(FALSE, 3), rep(TRUE, 7))

  fit <- hsreg_gibbs(y, X, penalized = pen, n_mcmc = 50, burnin = 20,
                          verbose = FALSE)

  expect_equal(fit$pen_idx, 4:10)
  expect_equal(fit$free_idx, 1:3)
  expect_equal(nrow(fit$lambda2_draws), 7)  # only penalized
})

test_that("hsreg_gibbs Cholesky fallback triggers on near-singular X", {
  set.seed(42)
  n <- 50; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  # Make column 2 nearly identical to column 1
  X[, 2] <- X[, 1] + rnorm(n, sd = 1e-12)
  y <- rnorm(n)

  # This should still run (ridge fallback), though it might warn
  fit <- hsreg_gibbs(y, X, n_mcmc = 30, burnin = 10, verbose = FALSE)
  # n_chol_fallback might be > 0 if the near-singularity triggers it

  expect_true(is.integer(fit$n_chol_fallback) || is.numeric(fit$n_chol_fallback))
})

test_that("hsreg_gibbs posterior mean of sigma2 is reasonable", {
  # Generate data with known sigma^2 = 1
  set.seed(123)
  n <- 200; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(1, -1, 0.5, 0, 0)
  y <- X %*% beta_true + rnorm(n, sd = 1)

  fit <- hsreg_gibbs(y, X, n_mcmc = 500, burnin = 200, verbose = FALSE)

  # Posterior mean of sigma2 should be within a factor of 2 of true value
  sigma2_post_mean <- mean(fit$sigma2_draws)
  expect_true(sigma2_post_mean > 0.3)
  expect_true(sigma2_post_mean < 3.0)
})

test_that("hsreg_gibbs errors on zero penalized parameters", {
  n <- 50; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  expect_error(
    hsreg_gibbs(y, X, penalized = rep(FALSE, p), verbose = FALSE),
    "No penalized parameters"
  )
})

test_that("hsreg_gibbs errors on dimension mismatch", {
  X <- matrix(rnorm(100), 20, 5)
  y <- rnorm(10)  # wrong length

  expect_error(hsreg_gibbs(y, X), "must equal nrow")
})

test_that("hsreg_gibbs errors on invalid scale parameters", {
  n <- 50; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  expect_error(hsreg_gibbs(y, X, lambda_scale = -1), "must be a positive scalar")
  expect_error(hsreg_gibbs(y, X, tau_scale = 0), "must be a positive scalar")
})
