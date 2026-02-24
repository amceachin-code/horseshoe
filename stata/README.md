# horseshoe (Stata)

Bayesian Linear Regression with the Horseshoe Prior -- Stata e-class command

## Overview

**horseshoe** is a Stata e-class estimation command implementing a Gibbs sampler
for Bayesian linear regression with the horseshoe prior (Carvalho, Polson &
Scott, 2010). It uses the auxiliary-variable sampler from Makalic & Schmidt
(2015, arXiv:1508.03884) with extensions for:

- **Selective penalization** -- per-parameter control via `penalized()` matrix
- **Convergence diagnostics** -- effective sample size (ESS) for every parameter
- **Numerical hardening** -- Cholesky fallback with ridge adjustment
- **Full postestimation support** -- `predict`, `test`, `lincom`, `margins`
- **CATE extraction** -- companion `horseshoe_cates` command
- **No dependencies** -- uses only base Stata and Mata

## Installation

Copy all five files from `ado/` to your personal ado directory (run `sysdir`
in Stata to find it), or point Stata at this folder:

```stata
adopath + "/path/to/hsreg/stata/ado"
```

Verify:

```stata
which horseshoe
help horseshoe
```

## Files

| File | Description |
|------|-------------|
| `ado/horseshoe.ado` | Main estimation command (Gibbs sampler + display) |
| `ado/horseshoe_p.ado` | Predict postestimation command |
| `ado/horseshoe_cates.ado` | CATE extraction command + Mata functions |
| `ado/_hs_utils.ado` | Shared Mata utilities (ESS, quantile pair) |
| `ado/horseshoe.sthlp` | Stata help file (`help horseshoe`) |

## Quick Example

```stata
sysuse auto, clear
horseshoe price mpg weight length turn, nmcmc(2000) burnin(500) seed(42)
matrix list e(b)
predict yhat, xb
```

### Selective Penalization

```stata
matrix pen = J(1, 20, 1)
forvalues j = 1/4 {
    matrix pen[1, `j'] = 0
}
horseshoe y treat1-treat4 x1-x16, penalized(pen) nmcmc(2000) seed(42)
```

### CATE Extraction

```stata
horseshoe y x1-x10 d1-d4 xd1-xd40, penalized(pen) saving(draws.dta) replace
horseshoe_cates, draws(draws.dta) xtest(X_test) nx(10) nd(4)
matrix list r(cate_hat)
```

## References

- Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe
  estimator for sparse signals. *Biometrika*, 97(2), 465-480.
- Makalic, E. & Schmidt, D. F. (2015). A simple sampler for the horseshoe
  estimator. *arXiv:1508.03884v4*.
