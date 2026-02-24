# hsreg

Bayesian Linear Regression with the Horseshoe Prior -- R package

## Overview

**hsreg** is an R package implementing a Gibbs sampler for Bayesian
linear regression with the horseshoe prior (Carvalho, Polson & Scott, 2010).
It uses the auxiliary-variable sampler from Makalic & Schmidt (2015,
arXiv:1508.03884) with extensions for:

- **Selective penalization** -- user-specified coefficients receive the
  horseshoe prior (shrinkage) while others receive a flat prior (no
  shrinkage). Leave treatment main effects unpenalized while shrinking
  high-dimensional interactions.
- **Convergence diagnostics** -- effective sample size (ESS) for every
  parameter via the initial positive sequence estimator (Geyer 1992).
- **Numerical hardening** -- Cholesky fallback with ridge adjustment for
  near-singular designs, NaN/Inf divergence detection with informative
  error messages.
- **CATE extraction** -- posterior CATEs with credible intervals for
  multi-arm treatment designs with treatment-covariate interactions.
- **Full S3 support** -- `print`, `summary`, `coef`, `vcov`, `predict`,
  `confint` methods that work like `lm`/`glm`.
- **Zero external dependencies** -- depends only on base R (`stats`).

## Installation

```r
# Install from GitHub
remotes::install_github("amceachin-code/horseshoe", subdir = "R")
```

## Quick Example

```r
library(hsreg)

# Sparse regression: 5 true signals in 50 predictors
set.seed(42)
n <- 200; p <- 50
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(rep(3, 5), rep(0, 45))
y <- X %*% beta_true + rnorm(n)

# Fit with horseshoe prior
fit <- hsreg(y, X, n_mcmc = 1000, burnin = 500, seed = 42)
print(fit)

# Predictions and credible intervals
predict(fit, type = "xb")[1:5]
confint(fit)[1:5, ]
```

### Selective Penalization

```r
# Treatment main effects (cols 11-12) unpenalized, everything else penalized
penalized <- c(rep(TRUE, 10), FALSE, FALSE, rep(TRUE, 38))
fit <- hsreg(y, X, penalized = penalized, n_mcmc = 1000, burnin = 500)
```

### CATE Extraction

```r
# Extract individual CATEs for multi-arm treatment design
# Design layout: [X(n_x), D(n_d), X:D(n_d * n_x)]
cates <- hsreg_cates(fit, X_test, n_x = 10, n_d = 4)
head(cates$cate_hat)  # posterior mean CATEs
head(cates$cate_lo)   # lower 95% CI
head(cates$cate_hi)   # upper 95% CI
```

## Public API

| Function | Description |
|----------|-------------|
| `hsreg()` | Fit Bayesian linear regression with horseshoe prior |
| `hsreg_gibbs()` | Raw Gibbs sampler (advanced) |
| `hsreg_cates()` | Extract CATEs from fitted object |
| `hsreg_cates_from_file()` | Extract CATEs from saved draws file |
| `hs_ess()` | Effective sample size for MCMC chain |
| `hs_quantile_pair()` | Paired quantile (single sort) |
| `print.hsreg()` | Display results table |
| `summary.hsreg()` | Summarize results |
| `coef.hsreg()` | Extract coefficient vector |
| `vcov.hsreg()` | Extract posterior VCV |
| `predict.hsreg()` | Predictions (xb, residuals, stdp) |
| `confint.hsreg()` | Credible intervals |

## References

- Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe
  estimator for sparse signals. *Biometrika*, 97(2), 465-480.
- Makalic, E. & Schmidt, D. F. (2015). A simple sampler for the horseshoe
  estimator. *arXiv:1508.03884v4*.
- Geyer, C. J. (1992). Practical Markov chain Monte Carlo. *Statistical
  Science*, 7(4), 473-483.

## License

GPL (>= 3)
