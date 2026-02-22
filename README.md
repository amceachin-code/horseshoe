# horseshoe

Bayesian linear regression with the horseshoe prior for Stata.

![Version](https://img.shields.io/badge/version-1.1.0-blue)
![Stata](https://img.shields.io/badge/Stata-16%2B-brightgreen)
![License](https://img.shields.io/badge/license-MIT-green)

## Overview

`horseshoe` fits a Bayesian linear regression model using the horseshoe prior ([Carvalho, Polson, and Scott 2010](https://doi.org/10.1093/biomet/asq017)), estimated via the auxiliary-variable Gibbs sampler of [Makalic and Schmidt (2015)](https://arxiv.org/abs/1508.03884).

The horseshoe prior is a continuous shrinkage prior that adapts to sparsity: it aggressively shrinks noise coefficients toward zero while allowing true signals to remain large. The key differentiator of this implementation is **per-parameter penalization control** via the `penalized()` option, which lets you specify exactly which coefficients receive the horseshoe prior and which receive a flat (improper) prior.

**The model:**

```
y = X * beta + epsilon,   epsilon ~ N(0, sigma^2)

For penalized coefficients (penalized[j] = 1):
    beta_j | lambda_j, tau, sigma^2  ~  N(0, lambda_j^2 * tau^2 * sigma^2)
    lambda_j  ~  C+(0, lambda_scale)
    tau       ~  C+(0, tau_scale)

For unpenalized coefficients (penalized[j] = 0):
    beta_j  ~  flat (improper) prior
```

The half-Cauchy priors are represented via scale-mixture of inverse-gamma distributions using auxiliary variables, yielding fully conjugate conditional posteriors for efficient Gibbs sampling.

## Features

- **Per-parameter penalization** via a `penalized()` matrix of 0/1 flags -- leave treatment main effects unpenalized while shrinking high-dimensional interactions
- **Convergence diagnostics** -- effective sample size (ESS) for every parameter using the initial positive sequence estimator ([Geyer 1992](https://doi.org/10.1214/ss/1177011137))
- **Reproducibility** -- `seed()` option sets the RNG state; `e(rng_state)` stores the full state for bitwise-identical replication
- **Full posterior access** -- `saving()` writes all MCMC draws to a `.dta` file with named beta columns for custom analysis
- **Standard Stata postestimation** -- posts `e(b)` and `e(V)`, enabling `predict`, `test`, `lincom`, `nlcom`, and `margins`
- **CATE extraction** -- companion command `horseshoe_cates` extracts individual-level conditional average treatment effects for multi-arm treatment designs
- **Numerical robustness** -- NaN/Inf detection at every iteration, Cholesky fallback with ridge adjustment for near-singular designs
- **Named betas in saved draws** -- saved `.dta` files use `b_varname` columns for readability
- **Replay support** -- running `horseshoe` without arguments re-displays the last results table

## Installation

Copy all 5 package files to any directory on your Stata adopath:

```
horseshoe.ado
horseshoe_p.ado
horseshoe_cates.ado
_hs_utils.ado
horseshoe.sthlp
```

Or point your adopath to the directory containing them:

```stata
adopath + "/path/to/horseshoe/files"
```

**Requirements:** Stata 16.0 or later. No external dependencies.

To verify installation:

```stata
which horseshoe
help horseshoe
```

## Quick Start

```stata
sysuse auto, clear
horseshoe price mpg weight length turn, nmcmc(2000) burnin(500) seed(42)
predict yhat, xb
```

## Usage Examples

### Standard horseshoe (all coefficients penalized)

```stata
sysuse auto, clear
horseshoe price mpg weight length turn, nmcmc(2000) burnin(500) verbose
```

### Selective penalization (treatment dummies unpenalized)

```stata
* Create penalization flag: 0 for treatment (cols 1-4), 1 for interactions (cols 5-20)
matrix pen = J(1, 20, 1)
matrix pen[1, 1] = 0
matrix pen[1, 2] = 0
matrix pen[1, 3] = 0
matrix pen[1, 4] = 0
horseshoe y treat1 treat2 treat3 treat4 x1-x16, penalized(pen)
```

### Reproducible results with seed()

```stata
horseshoe price mpg weight length, seed(12345)
matrix list e(b)

* Run again with same seed -- identical results
horseshoe price mpg weight length, seed(12345)

* Or restore from saved RNG state
local saved_state = e(rng_state)
horseshoe price mpg weight length, seed(`saved_state')
```

### Saving and analyzing posterior draws

```stata
horseshoe price mpg weight length, saving(posterior_draws.dta) replace seed(42)

* Load draws and compute custom summaries
preserve
use posterior_draws.dta, clear
summarize sigma2 tau2

* Posterior probability that b_mpg > 0
count if b_mpg > 0
display "Pr(beta_mpg > 0) = " r(N) / _N
restore
```

### Checking convergence diagnostics

```stata
horseshoe price mpg weight length
display "Min ESS = " e(ess_min)
display "ESS sigma^2 = " e(ess_sigma2)
display "ESS tau^2 = " e(ess_tau2)
matrix list e(ess_beta)
```

### Predict and postestimation

```stata
sysuse auto, clear
horseshoe price mpg weight length, nmcmc(1000) seed(42)

* Fitted values
predict yhat, xb
* Residuals
predict resid, residuals
* Standard errors of prediction
predict se_yhat, stdp

* Test a linear hypothesis
test mpg = weight
* Linear combination
lincom mpg + weight

* Replay results
horseshoe
```

### CATE extraction for multi-arm treatment designs

```stata
* After fitting horseshoe with saving()
horseshoe_cates, draws(posterior.dta) xtest(X_test) xsd(X_sd) nx(54) nd(4)

* Results returned as r-class matrices
matrix list r(cate_hat)   // n_test x n_d posterior mean iCATEs
matrix list r(cate_lo)    // n_test x n_d lower CI bounds
matrix list r(cate_hi)    // n_test x n_d upper CI bounds
```

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `penalized(matname)` | matrix name | all 1s | 1 x p matrix of 0/1 flags: 1 = horseshoe prior, 0 = flat prior |
| `nmcmc(#)` | integer | 1000 | Number of post-burnin MCMC draws to store |
| `burnin(#)` | integer | 500 | Number of burnin iterations to discard |
| `thin(#)` | integer | 1 | Thinning interval (keep every thin-th draw) |
| `lambdascale(#)` | real | 1 | Scale of the half-Cauchy prior on local shrinkage lambda_j |
| `tauscale(#)` | real | 1 | Scale of the half-Cauchy prior on global shrinkage tau |
| `saving(filename)` | string | -- | Save full posterior draws to a .dta file |
| `replace` | flag | off | Allow overwriting an existing `saving()` file |
| `seed(#\|string)` | integer or string | -- | RNG seed (integer) or full RNG state string |
| `level(#)` | real | 95 | Credible interval level (percentage) |
| `verbose` | flag | off | Print progress every 500 iterations |

## Stored Results

### Scalars

| Name | Description |
|---|---|
| `e(N)` | Number of observations |
| `e(p)` | Number of predictors |
| `e(p_pen)` | Number of penalized predictors |
| `e(nmcmc)` | Number of stored MCMC draws |
| `e(sigma2)` | Posterior mean of sigma^2 |
| `e(tau2)` | Posterior mean of tau^2 |
| `e(lambda_scale)` | Lambda scale used |
| `e(tau_scale)` | Tau scale used |
| `e(burnin)` | Number of burnin iterations |
| `e(thin)` | Thinning interval |
| `e(n_chol_fallback)` | Iterations requiring Cholesky ridge adjustment |
| `e(ess_min)` | Minimum ESS across all parameters |
| `e(ess_sigma2)` | ESS for sigma^2 |
| `e(ess_tau2)` | ESS for tau^2 |

### Macros

| Name | Description |
|---|---|
| `e(cmd)` | `horseshoe` |
| `e(cmdline)` | Full command as typed |
| `e(predict)` | `horseshoe_p` |
| `e(depvar)` | Name of dependent variable |
| `e(seed)` | Seed specified (empty if not specified) |
| `e(rng_state)` | Full RNG state captured before sampling |

### Matrices

| Name | Dimensions | Description |
|---|---|---|
| `e(b)` | 1 x p | Posterior means of beta |
| `e(V)` | p x p | Posterior variance-covariance matrix of beta |
| `e(b_sd)` | 1 x p | Posterior standard deviations of beta |
| `e(b_lower)` | 1 x p | Lower credible interval bounds |
| `e(b_upper)` | 1 x p | Upper credible interval bounds |
| `e(penalized)` | 1 x p | Penalization flags (0/1) |
| `e(ess_beta)` | 1 x p | Effective sample sizes for each beta_j |

## Convergence Assessment

`horseshoe` automatically computes the effective sample size (ESS) for all sampled parameters using the initial positive sequence estimator (Geyer 1992).

**ESS guidelines:**

- **ESS > 400**: Generally adequate for reliable posterior summaries (means, SDs, 95% CIs)
- **ESS 100--400**: Sufficient for posterior means; credible interval tails may be imprecise
- **ESS < 100**: Poor mixing -- increase `nmcmc()` or `burnin()`, or check for near-collinearity

**tau^2 is typically the mixing bottleneck.** The global shrinkage parameter mixes slowest. If `e(ess_tau2)` is much lower than beta ESS values, the chain needs more iterations.

**Tuning guidance:**

- Defaults (`nmcmc(1000) burnin(500)`) are adequate for moderate p (up to ~100 predictors)
- For p > 200, try `nmcmc(2000) burnin(1000)`
- Thinning (`thin()`) is usually unnecessary -- it does not improve ESS per unit of compute time

A warning is displayed automatically when the minimum ESS falls below 100.

## CATE Extraction

The companion command `horseshoe_cates` extracts individual-level conditional average treatment effects (iCATEs) from posterior draws for multi-arm treatment designs.

**Stata command:**

```stata
horseshoe_cates, draws(filename) xtest(matname) xsd(matname) nx(#) nd(#) [level(#)]
```

**Mata function** (for use in pipeline scripts):

```stata
run "horseshoe_cates.ado"
mata: _hs_extract_cates(beta_all, X_test, X_sd, n_x, n_d, lo_prob, hi_prob)
```

**Column layout assumption** (must match the design matrix used during estimation):

| Columns | Content |
|---|---|
| 1 to n_x | X covariate main effects |
| n_x+1 to n_x+n_d | D treatment dummies (e.g., D_b, D_c, D_d, D_e) |
| n_x+n_d+1 to end | X:D interactions, blocked by arm (n_x columns per arm) |

## Package Files

| File | Lines | Description |
|---|---|---|
| `horseshoe.ado` | 1,309 | Main estimation command: syntax parsing, Mata Gibbs sampler, e-class posting |
| `horseshoe_cates.ado` | 303 | CATE extraction for multi-arm designs (Stata command + Mata function) |
| `horseshoe_p.ado` | 70 | Predict postestimation: xb, residuals, stdp |
| `_hs_utils.ado` | 86 | Shared Mata utilities (_hs_quantile_pair) |
| `horseshoe.sthlp` | 502 | Stata help file |
| **Total** | **2,270** | |

## References

- Carvalho, C. M., N. G. Polson, and J. G. Scott (2010). "The horseshoe estimator for sparse signals." *Biometrika* 97(2): 465--480. [doi:10.1093/biomet/asq017](https://doi.org/10.1093/biomet/asq017)
- Makalic, E. and D. F. Schmidt (2015). "A simple sampler for the horseshoe estimator." [arXiv:1508.03884](https://arxiv.org/abs/1508.03884)
- Geyer, C. J. (1992). "Practical Markov chain Monte Carlo." *Statistical Science* 7(4): 473--483. [doi:10.1214/ss/1177011137](https://doi.org/10.1214/ss/1177011137)

## License

MIT

## Author

Andrew McEachin
