# ===========================================================================
# tests/testthat/test-utils.R — Tests for hs_ess() and hs_quantile_pair()
# ===========================================================================

# ---------------------------------------------------------------------------
# hs_ess()
# ---------------------------------------------------------------------------

test_that("hs_ess returns ESS in [1, n] range", {
  set.seed(42)
  chain <- cumsum(rnorm(500))  # random walk — highly autocorrelated
  ess <- hs_ess(chain)
  expect_true(ess >= 1)
  expect_true(ess <= 500)
})

test_that("hs_ess on constant chain returns n", {
  # A constant chain has zero variance — all draws are identical,

  # so ESS should be n (no information loss from autocorrelation)
  chain <- rep(3.14, 200)
  ess <- hs_ess(chain)
  expect_equal(ess, 200)
})

test_that("hs_ess on IID chain returns approximately n", {
  # IID draws should have ESS close to n
  set.seed(123)
  chain <- rnorm(1000)
  ess <- hs_ess(chain)
  # Should be within 30% of n for IID
  expect_true(ess > 700)
  expect_true(ess <= 1000)
})

test_that("hs_ess on highly autocorrelated chain returns << n", {
  # AR(1) with rho = 0.99 — very slow mixing
  set.seed(99)
  n <- 2000
  chain <- numeric(n)
  chain[1] <- rnorm(1)
  for (i in 2:n) chain[i] <- 0.99 * chain[i-1] + rnorm(1, sd = 0.1)
  ess <- hs_ess(chain)
  # ESS should be much less than n (< 10% of n for rho=0.99)
  expect_true(ess < n * 0.15)
})

test_that("hs_ess handles short chains gracefully", {
  expect_equal(hs_ess(c(1)), 1)
  expect_equal(hs_ess(c(1, 2)), 2)
  expect_equal(hs_ess(c(1, 2, 3)), 3)
})

test_that("hs_ess rejects non-numeric input", {
  expect_error(hs_ess("abc"), "'chain' must be a numeric")
})


# ---------------------------------------------------------------------------
# hs_quantile_pair()
# ---------------------------------------------------------------------------

test_that("hs_quantile_pair matches stats::quantile on random vectors", {
  set.seed(42)
  x <- rnorm(1000)

  # 95% CI
  result <- hs_quantile_pair(x, 0.025, 0.975)
  expected_lo <- unname(quantile(x, 0.025))
  expected_hi <- unname(quantile(x, 0.975))

  expect_equal(result[["lower"]], expected_lo, tolerance = 1e-10)
  expect_equal(result[["upper"]], expected_hi, tolerance = 1e-10)
})

test_that("hs_quantile_pair matches stats::quantile at various levels", {
  set.seed(7)
  x <- runif(500)

  # 90% CI
  result <- hs_quantile_pair(x, 0.05, 0.95)
  expect_equal(result[["lower"]], unname(quantile(x, 0.05)), tolerance = 1e-10)
  expect_equal(result[["upper"]], unname(quantile(x, 0.95)), tolerance = 1e-10)

  # Median as both bounds (degenerate case)
  result2 <- hs_quantile_pair(x, 0.5, 0.5)
  expect_equal(result2[["lower"]], unname(quantile(x, 0.5)), tolerance = 1e-10)
  expect_equal(result2[["upper"]], unname(quantile(x, 0.5)), tolerance = 1e-10)
})

test_that("hs_quantile_pair returns named vector", {
  result <- hs_quantile_pair(1:10, 0.1, 0.9)
  expect_named(result, c("lower", "upper"))
})

test_that("hs_quantile_pair handles single element", {
  result <- hs_quantile_pair(42, 0.025, 0.975)
  expect_equal(result[["lower"]], 42)
  expect_equal(result[["upper"]], 42)
})

test_that("hs_quantile_pair rejects empty vector", {
  expect_error(hs_quantile_pair(numeric(0), 0.025, 0.975))
})

test_that("hs_quantile_pair rejects non-numeric input", {
  expect_error(hs_quantile_pair("abc", 0.025, 0.975), "'v' must be a numeric")
})
