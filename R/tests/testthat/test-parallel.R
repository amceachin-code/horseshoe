# ===========================================================================
# tests/testthat/test-parallel.R — Tests for multi-chain parallel support
# ===========================================================================

# --- Helper: small reproducible data used across tests ---
make_test_data <- function(seed = 42, n = 100, p = 10) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", 1:p)
  beta_true <- c(rep(2, 3), rep(0, p - 3))
  y <- X %*% beta_true + rnorm(n)
  list(y = y, X = X)
}


# ---------------------------------------------------------------------------
# Backward compatibility: chains=1 is bit-for-bit identical to default
# ---------------------------------------------------------------------------

test_that("chains=1 produces identical results to default (no chains arg)", {
  dat <- make_test_data()

  # Default call (original code path)
  fit_default <- hsreg(dat$y, dat$X, n_mcmc = 100, burnin = 50,
                       verbose = FALSE, seed = 42)

  # Explicit chains=1
  fit_single <- hsreg(dat$y, dat$X, n_mcmc = 100, burnin = 50,
                      verbose = FALSE, seed = 42, chains = 1)

  expect_identical(fit_default$beta_draws, fit_single$beta_draws)
  expect_identical(fit_default$sigma2_draws, fit_single$sigma2_draws)
  expect_identical(fit_default$tau2_draws, fit_single$tau2_draws)
  expect_identical(fit_default$coefficients, fit_single$coefficients)
})


# ---------------------------------------------------------------------------
# Multi-chain reproducibility: same seed + chains = same draws
# ---------------------------------------------------------------------------

test_that("multi-chain is reproducible with same seed", {
  dat <- make_test_data()

  fit1 <- hsreg(dat$y, dat$X, n_mcmc = 50, burnin = 25,
                verbose = FALSE, seed = 123, chains = 2, cores = 1)
  fit2 <- hsreg(dat$y, dat$X, n_mcmc = 50, burnin = 25,
                verbose = FALSE, seed = 123, chains = 2, cores = 1)

  expect_identical(fit1$beta_draws, fit2$beta_draws)
  expect_identical(fit1$sigma2_draws, fit2$sigma2_draws)
  expect_identical(fit1$tau2_draws, fit2$tau2_draws)
})


# ---------------------------------------------------------------------------
# Draw dimensions: total_draws == chains * n_mcmc
# ---------------------------------------------------------------------------

test_that("draw dimensions are correct for multi-chain", {
  dat <- make_test_data()
  n_mcmc <- 50
  chains <- 3

  fit <- hsreg(dat$y, dat$X, n_mcmc = n_mcmc, burnin = 25,
               verbose = FALSE, seed = 42, chains = chains, cores = 1)

  expect_equal(ncol(fit$beta_draws), chains * n_mcmc)
  expect_equal(length(fit$sigma2_draws), chains * n_mcmc)
  expect_equal(length(fit$tau2_draws), chains * n_mcmc)
  expect_equal(ncol(fit$lambda2_draws), chains * n_mcmc)
  expect_equal(fit$total_draws, chains * n_mcmc)
  expect_equal(fit$chains, chains)
})


# ---------------------------------------------------------------------------
# R-hat computed for multi-chain, NULL for single chain
# ---------------------------------------------------------------------------

test_that("R-hat is computed for multi-chain and NULL for single chain", {
  dat <- make_test_data()

  # Single chain: R-hat should be NULL
  fit_single <- hsreg(dat$y, dat$X, n_mcmc = 50, burnin = 25,
                      verbose = FALSE, seed = 42, chains = 1)
  expect_null(fit_single$rhat_beta)
  expect_null(fit_single$rhat_sigma2)
  expect_null(fit_single$rhat_tau2)
  expect_null(fit_single$rhat_max)

  # Multi-chain: R-hat should be non-NULL, numeric, >= 1.0
  fit_multi <- hsreg(dat$y, dat$X, n_mcmc = 50, burnin = 25,
                     verbose = FALSE, seed = 42, chains = 2, cores = 1)
  expect_true(!is.null(fit_multi$rhat_beta))
  expect_true(is.numeric(fit_multi$rhat_beta))
  expect_length(fit_multi$rhat_beta, ncol(dat$X))
  expect_true(all(fit_multi$rhat_beta >= 0.9))  # should be near 1.0

  expect_true(is.numeric(fit_multi$rhat_sigma2))
  expect_true(is.numeric(fit_multi$rhat_tau2))
  expect_true(is.numeric(fit_multi$rhat_max))
  expect_true(fit_multi$rhat_max >= 0.9)
})


# ---------------------------------------------------------------------------
# R-hat unit tests: converged chains → ~1.0, divergent → >>1
# ---------------------------------------------------------------------------

test_that(".hs_rhat returns ~1.0 for converged chains", {
  set.seed(42)
  # Two chains from the same normal distribution → R-hat ≈ 1.0
  chain1 <- rnorm(1000, mean = 5, sd = 1)
  chain2 <- rnorm(1000, mean = 5, sd = 1)
  rhat <- .hs_rhat(list(chain1, chain2))
  expect_true(rhat < 1.05)
  expect_true(rhat >= 1.0 || abs(rhat - 1.0) < 0.01)
})

test_that(".hs_rhat returns >>1 for divergent chains", {
  set.seed(42)
  # Two chains with very different means → R-hat >> 1
  chain1 <- rnorm(500, mean = 0, sd = 1)
  chain2 <- rnorm(500, mean = 100, sd = 1)
  rhat <- .hs_rhat(list(chain1, chain2))
  expect_true(rhat > 1.5)
})

test_that(".hs_rhat_vector works for a matrix of parameters", {
  set.seed(42)
  p <- 5
  n_mcmc <- 200

  # Two chain matrices, both from same distribution
  mat1 <- matrix(rnorm(p * n_mcmc), p, n_mcmc)
  mat2 <- matrix(rnorm(p * n_mcmc), p, n_mcmc)
  rhat_vec <- .hs_rhat_vector(list(mat1, mat2))

  expect_length(rhat_vec, p)
  expect_true(all(rhat_vec < 1.1))
})


# ---------------------------------------------------------------------------
# S3 methods work with multi-chain fits
# ---------------------------------------------------------------------------

test_that("S3 methods work with multi-chain fit", {
  dat <- make_test_data()
  fit <- hsreg(dat$y, dat$X, n_mcmc = 50, burnin = 25,
               verbose = FALSE, seed = 42, chains = 2, cores = 1)

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
  expect_equal(s$chains, 2L)
  expect_equal(s$total_draws, 100L)
  expect_true(!is.null(s$rhat_max))
  expect_no_error(capture.output(print(s)))
})


# ---------------------------------------------------------------------------
# hsreg_cates works with multi-chain fit
# ---------------------------------------------------------------------------

test_that("hsreg_cates works with multi-chain fit", {
  set.seed(42)
  n <- 100; n_x <- 5; n_d <- 2
  p <- n_x + n_d + n_d * n_x  # 5 + 2 + 10 = 17

  X_design <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  X_test <- matrix(rnorm(20 * n_x), 20, n_x)

  # Fit with 2 chains
  fit <- hsreg(y, X_design, n_mcmc = 50, burnin = 25,
               verbose = FALSE, seed = 42, chains = 2, cores = 1)

  # n_draws should be total_draws = 2 * 50 = 100
  cates <- hsreg_cates(fit, X_test = X_test, n_x = n_x, n_d = n_d)
  expect_equal(cates$n_draws, 100)
  expect_equal(dim(cates$cate_hat), c(20, n_d))
  expect_equal(dim(cates$cate_lo), c(20, n_d))
  expect_equal(dim(cates$cate_hi), c(20, n_d))
})


# ---------------------------------------------------------------------------
# Parallel == sequential draws: same seed/chains, different cores → identical
# ---------------------------------------------------------------------------

test_that("parallel and sequential produce identical draws", {
  dat <- make_test_data()

  # Sequential (cores=1)
  fit_seq <- hsreg(dat$y, dat$X, n_mcmc = 50, burnin = 25,
                   verbose = FALSE, seed = 42, chains = 2, cores = 1)

  # Parallel (cores=2)
  fit_par <- hsreg(dat$y, dat$X, n_mcmc = 50, burnin = 25,
                   verbose = FALSE, seed = 42, chains = 2, cores = 2)

  expect_identical(fit_seq$beta_draws, fit_par$beta_draws)
  expect_identical(fit_seq$sigma2_draws, fit_par$sigma2_draws)
  expect_identical(fit_seq$tau2_draws, fit_par$tau2_draws)
})


# ---------------------------------------------------------------------------
# Saving + multi-chain round-trip
# ---------------------------------------------------------------------------

test_that("saving works correctly with multi-chain", {
  dat <- make_test_data(n = 50, p = 5)
  tmpfile <- tempfile(fileext = ".rds")
  on.exit(unlink(tmpfile))

  n_mcmc <- 50
  chains <- 2

  fit <- hsreg(dat$y, dat$X, n_mcmc = n_mcmc, burnin = 25,
               saving = tmpfile, verbose = FALSE, seed = 42,
               chains = chains, cores = 1)

  # File should exist
  expect_true(file.exists(tmpfile))

  # Load and check dimensions
  draws <- readRDS(tmpfile)
  expect_true(is.data.frame(draws))
  expect_equal(nrow(draws), chains * n_mcmc)  # total_draws rows
  expect_true("draw" %in% colnames(draws))
  expect_true("sigma2" %in% colnames(draws))
  expect_true("tau2" %in% colnames(draws))
})

test_that("hsreg_cates_from_file works with multi-chain saved draws", {
  set.seed(42)
  n <- 80; n_x <- 4; n_d <- 2
  p <- n_x + n_d + n_d * n_x  # 4 + 2 + 8 = 14

  X_design <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  X_test <- matrix(rnorm(10 * n_x), 10, n_x)

  tmpfile <- tempfile(fileext = ".rds")
  on.exit(unlink(tmpfile))

  # Fit and save with 2 chains
  fit <- hsreg(y, X_design, n_mcmc = 30, burnin = 15,
               saving = tmpfile, verbose = FALSE, seed = 42,
               chains = 2, cores = 1)

  # Load CATEs from file — should handle 60 total draws

  cates <- hsreg_cates_from_file(tmpfile, X_test = X_test,
                                 n_x = n_x, n_d = n_d)
  expect_equal(cates$n_draws, 60)  # 2 chains * 30 draws
  expect_equal(dim(cates$cate_hat), c(10, n_d))
})


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("invalid chains/cores values produce errors", {
  dat <- make_test_data(n = 50, p = 5)

  expect_error(
    hsreg(dat$y, dat$X, n_mcmc = 10, burnin = 5, verbose = FALSE, chains = 0),
    "'chains' must be a positive integer"
  )
  expect_error(
    hsreg(dat$y, dat$X, n_mcmc = 10, burnin = 5, verbose = FALSE, chains = -1),
    "'chains' must be a positive integer"
  )
  expect_error(
    hsreg(dat$y, dat$X, n_mcmc = 10, burnin = 5, verbose = FALSE, cores = -1),
    "'cores' must be a positive integer"
  )
  expect_error(
    hsreg(dat$y, dat$X, n_mcmc = 10, burnin = 5, verbose = FALSE, cores = 0),
    "'cores' must be a positive integer"
  )
})
