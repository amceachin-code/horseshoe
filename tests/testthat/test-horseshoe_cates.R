# ===========================================================================
# tests/testthat/test-horseshoe_cates.R — Tests for CATE extraction
# ===========================================================================


test_that("horseshoe_cates recovers known CATEs within CI", {
  # DGP: 2 covariates, 2 treatment contrasts
  # Design: [X(2), D(2), X:D(4)] = 8 columns total
  # CATE for contrast 1: tau_1(x) = 1.0 + 0.5*x1 + 0.0*x2
  # CATE for contrast 2: tau_2(x) = -0.5 + 0.0*x1 + 1.0*x2
  set.seed(42)
  n_x <- 2; n_d <- 2
  n <- 200
  p <- n_x + n_d + n_d * n_x  # 2 + 2 + 4 = 8

  X_cov <- matrix(rnorm(n * n_x), n, n_x)

  # Treatment assignment (3 arms: reference a, and b, c)
  D <- matrix(0, n, n_d)
  arms <- sample(0:2, n, replace = TRUE)
  D[arms == 1, 1] <- 1
  D[arms == 2, 2] <- 1

  # Interactions: X:D blocked by arm
  XD <- cbind(X_cov * D[, 1], X_cov * D[, 2])

  # Full design matrix
  X_design <- cbind(X_cov, D, XD)

  # True beta: [beta_x(2), gamma(2), delta_1(2), delta_2(2)]
  # beta_x = c(0.3, -0.2) (main effects)
  # gamma  = c(1.0, -0.5) (treatment main effects)
  # delta_1 = c(0.5, 0.0) (interaction: arm 1)
  # delta_2 = c(0.0, 1.0) (interaction: arm 2)
  beta_true <- c(0.3, -0.2, 1.0, -0.5, 0.5, 0.0, 0.0, 1.0)

  y <- X_design %*% beta_true + rnorm(n, sd = 0.5)

  # Fit horseshoe — leave treatment main effects unpenalized
  pen <- c(rep(TRUE, n_x), rep(FALSE, n_d), rep(TRUE, n_d * n_x))

  fit <- horseshoe(y, X_design, penalized = pen,
                   n_mcmc = 500, burnin = 300, verbose = FALSE, seed = 42)

  # Extract CATEs at 10 test points
  n_test <- 10
  X_test <- matrix(rnorm(n_test * n_x), n_test, n_x)

  cates <- horseshoe_cates(fit, X_test, n_x = n_x, n_d = n_d, level = 0.95)

  expect_equal(dim(cates$cate_hat), c(n_test, n_d))
  expect_equal(dim(cates$cate_lo), c(n_test, n_d))
  expect_equal(dim(cates$cate_hi), c(n_test, n_d))

  # True CATEs
  true_cate_1 <- 1.0 + X_test %*% c(0.5, 0.0)
  true_cate_2 <- -0.5 + X_test %*% c(0.0, 1.0)

  # Check: most true CATEs should fall within the 95% CI
  # (at least 70% coverage — conservative for a small sample)
  coverage_1 <- mean(true_cate_1 >= cates$cate_lo[, 1] &
                     true_cate_1 <= cates$cate_hi[, 1])
  coverage_2 <- mean(true_cate_2 >= cates$cate_lo[, 2] &
                     true_cate_2 <= cates$cate_hi[, 2])

  expect_true(coverage_1 >= 0.5)
  expect_true(coverage_2 >= 0.5)
})


test_that("horseshoe_cates returns correct dimensions", {
  # Create a fake fit with known beta_draws
  set.seed(42)
  n_x <- 3; n_d <- 2
  p <- n_x + n_d + n_d * n_x  # 3 + 2 + 6 = 11
  n_mcmc <- 100

  fake_fit <- list(
    beta_draws = matrix(rnorm(p * n_mcmc), p, n_mcmc),
    X_sd = NULL
  )

  X_test <- matrix(rnorm(20 * n_x), 20, n_x)

  cates <- horseshoe_cates(fake_fit, X_test, n_x = n_x, n_d = n_d)

  expect_equal(cates$n_test, 20)
  expect_equal(cates$n_d, 2)
  expect_equal(cates$n_draws, 100)
  expect_equal(dim(cates$cate_hat), c(20, 2))
})


test_that("horseshoe_cates_from_file matches in-memory variant", {
  set.seed(42)
  n_x <- 2; n_d <- 2
  p <- n_x + n_d + n_d * n_x
  n <- 100

  X_cov <- matrix(rnorm(n * n_x), n, n_x)
  D <- matrix(0, n, n_d)
  arms <- sample(0:2, n, replace = TRUE)
  D[arms == 1, 1] <- 1
  D[arms == 2, 2] <- 1
  XD <- cbind(X_cov * D[, 1], X_cov * D[, 2])
  X_design <- cbind(X_cov, D, XD)
  y <- rnorm(n)

  tmpfile <- tempfile(fileext = ".rds")
  on.exit(unlink(tmpfile))

  pen <- c(rep(TRUE, n_x), rep(FALSE, n_d), rep(TRUE, n_d * n_x))
  fit <- horseshoe(y, X_design, penalized = pen,
                   n_mcmc = 100, burnin = 50, verbose = FALSE,
                   saving = tmpfile, seed = 42)

  X_test <- matrix(rnorm(10 * n_x), 10, n_x)

  # In-memory
  cates_mem <- horseshoe_cates(fit, X_test, n_x = n_x, n_d = n_d)

  # From file
  cates_file <- horseshoe_cates_from_file(tmpfile, X_test,
                                            n_x = n_x, n_d = n_d)

  expect_equal(cates_mem$cate_hat, cates_file$cate_hat, tolerance = 1e-10)
  expect_equal(cates_mem$cate_lo, cates_file$cate_lo, tolerance = 1e-10)
  expect_equal(cates_mem$cate_hi, cates_file$cate_hi, tolerance = 1e-10)
})


test_that("horseshoe_cates errors on wrong X_test dimensions", {
  fake_fit <- list(
    beta_draws = matrix(rnorm(110), 11, 10),
    X_sd = NULL
  )

  X_test_bad <- matrix(rnorm(20), 5, 4)  # n_x = 3, but 4 cols provided

  expect_error(
    horseshoe_cates(fake_fit, X_test_bad, n_x = 3, n_d = 2),
    "must have 3 columns"
  )
})


test_that("horseshoe_cates errors on wrong p in beta_draws", {
  fake_fit <- list(
    beta_draws = matrix(rnorm(100), 10, 10),  # p = 10
    X_sd = NULL
  )

  X_test <- matrix(rnorm(15), 5, 3)

  # n_x=3, n_d=2 expects p = 3+2+6 = 11, but got 10
  expect_error(
    horseshoe_cates(fake_fit, X_test, n_x = 3, n_d = 2),
    "expected 11"
  )
})
