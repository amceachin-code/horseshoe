# horseshoe

Bayesian Linear Regression with the Horseshoe Prior

Gibbs sampler for Bayesian linear regression with the horseshoe prior (Makalic
& Schmidt 2015). Features per-parameter penalization control, seed-based
reproducibility, ESS convergence diagnostics, Cholesky fallback and NaN/Inf
detection, and multi-arm CATE extraction. Available as both a Stata e-class
command and an R package.

## Implementations

| Language | Folder | Package Name | Install |
|----------|--------|-------------|---------|
| **R** | [`R/`](R/) | `hsreg` | `remotes::install_github("amceachin-code/horseshoe", subdir = "R")` |
| **Stata** | [`stata/`](stata/) | `horseshoe` | Copy `ado/` files to adopath |

Both implementations use the same algorithm -- the auxiliary-variable sampler
from Makalic & Schmidt (2015, arXiv:1508.03884) -- and produce statistically
equivalent results.

## Features

- **Selective penalization** -- horseshoe shrinkage on user-specified parameters,
  flat priors on the rest. Leave treatment main effects unpenalized while
  regularizing high-dimensional interactions.
- **Convergence diagnostics** -- effective sample size (ESS) via the initial
  positive sequence estimator (Geyer 1992) for every parameter.
- **Numerical hardening** -- Cholesky fallback with ridge adjustment for
  near-singular designs, NaN/Inf divergence detection.
- **CATE extraction** -- posterior conditional average treatment effects with
  credible intervals for multi-arm treatment designs.
- **Full postestimation support** -- R: `print`, `summary`, `coef`, `vcov`,
  `predict`, `confint`. Stata: `predict`, `test`, `lincom`, `margins`.
- **Zero external dependencies** -- R: base R + `stats`. Stata: base Stata + Mata.

## References

- Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe
  estimator for sparse signals. *Biometrika*, 97(2), 465-480.
- Makalic, E. & Schmidt, D. F. (2015). A simple sampler for the horseshoe
  estimator. *arXiv:1508.03884v4*.
- Geyer, C. J. (1992). Practical Markov chain Monte Carlo. *Statistical
  Science*, 7(4), 473-483.

## License

GPL (>= 3)
