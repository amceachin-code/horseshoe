*! horseshoe_cates.ado — CATE extraction from horseshoe posterior draws
*! Version 1.1.0 — 2026-02-21
*!
*! Provides both:
*!   1. Stata command: horseshoe_cates, draws(file) xtest(mat) xsd(mat) nx(#) nd(#) [level(#)]
*!   2. Mata function: _hs_extract_cates(beta_all, X_test, X_sd, n_x, n_d, lo, hi)
*!
*! Author: Andrew McEachin
*! License: MIT


* ===========================================================================
* STATA COMMAND WRAPPER
*
* horseshoe_cates, draws(filename) xtest(matname) xsd(matname)
*                  nx(#) nd(#) [level(#)]
*
* Loads posterior draws from a .dta file saved by horseshoe, saving(),
* extracts the beta columns, and calls _hs_extract_cates(). Results are
* returned as r-class matrices:
*   r(cate_hat)  — n_test × n_d posterior mean iCATEs
*   r(cate_lo)   — n_test × n_d lower credible interval bounds
*   r(cate_hi)   — n_test × n_d upper credible interval bounds
* ===========================================================================
program define horseshoe_cates, rclass
    version 16.0

    * Load shared Mata utilities (_hs_quantile_pair, etc.)
    capture findfile _hs_utils.ado
    if _rc {
        display as error "Cannot find _hs_utils.ado on the adopath."
        display as error ///
            "Ensure all horseshoe package files are in the same directory."
        exit 601
    }
    run "`r(fn)'"

    * Parse syntax
    syntax, DRaws(string) Xtest(name) Xsd(name) NX(integer) ND(integer) ///
            [Level(cilevel)]

    * -------------------------------------------------------------------
    * INPUT VALIDATION
    * -------------------------------------------------------------------

    * Check draw file exists
    confirm file "`draws'"

    * Check matrices exist
    confirm matrix `xtest'
    confirm matrix `xsd'

    * Validate nx and nd are positive
    if `nx' < 1 {
        display as error "nx() must be a positive integer, got `nx'"
        exit 198
    }
    if `nd' < 1 {
        display as error "nd() must be a positive integer, got `nd'"
        exit 198
    }

    * Validate X_test dimensions: should be n_test × nx
    if colsof(`xtest') != `nx' {
        display as error ///
            "xtest() must have `nx' columns (nx), got " colsof(`xtest')
        exit 503
    }

    * Expected total parameters: nx + nd + nd*nx
    local p_expected = `nx' + `nd' + `nd' * `nx'

    * Validate X_sd dimensions: should be 1 × p_expected
    if colsof(`xsd') != `p_expected' {
        display as error ///
            "xsd() must have `p_expected' columns (nx + nd + nd*nx), " ///
            "got " colsof(`xsd')
        exit 503
    }

    * Compute CI bounds from the level option
    local alpha = (100 - `level') / 100
    local lower_pct = `alpha' / 2
    local upper_pct = 1 - `alpha' / 2

    * -------------------------------------------------------------------
    * LOAD DRAWS AND EXTRACT BETA MATRIX
    * -------------------------------------------------------------------
    preserve
    quietly use "`draws'", clear

    * The draw file has columns: draw, sigma2, tau2, b_varname1, ..., b_varnameP
    * We need the beta columns (everything after draw, sigma2, tau2)
    quietly describe, short
    local nvars = r(k)
    local n_beta = `nvars' - 3

    if `n_beta' != `p_expected' {
        display as error ///
            "Draw file has `n_beta' beta columns but expected `p_expected' " ///
            "(nx + nd + nd*nx = `nx' + `nd' + " `nd' * `nx' ")"
        restore
        exit 503
    }

    * Get the variable names for the beta columns (columns 4 onward)
    quietly describe, varlist
    local allvars `r(varlist)'
    local betavars : list allvars - draw
    local betavars : list betavars - sigma2
    local betavars : list betavars - tau2

    * Convert beta columns to a Mata matrix
    mata: st_local("n_draws", strofreal(st_nobs()))
    mata: _hs_cates_from_draws("`betavars'", "`xtest'", "`xsd'", ///
        `nx', `nd', `lower_pct', `upper_pct')

    restore

    * -------------------------------------------------------------------
    * POST R-CLASS RESULTS
    * -------------------------------------------------------------------
    tempname cate_hat cate_lo cate_hi
    matrix `cate_hat' = __hs_cate_hat
    matrix `cate_lo'  = __hs_cate_lo
    matrix `cate_hi'  = __hs_cate_hi

    * Capture dimensions before return matrix moves the tempnames
    local n_test_val = rowsof(`cate_hat')

    return matrix cate_hat = `cate_hat'
    return matrix cate_lo  = `cate_lo'
    return matrix cate_hi  = `cate_hi'
    return scalar n_test   = `n_test_val'
    return scalar n_d      = `nd'
    return scalar n_draws  = `n_draws'
    return scalar level    = `level'

    * Clean up temporary Stata matrices
    capture matrix drop __hs_cate_hat __hs_cate_lo __hs_cate_hi

    display as text "horseshoe_cates: extracted iCATEs for " ///
        as result `n_test_val' as text " test observations, " ///
        as result `nd' as text " contrasts, " ///
        as result `n_draws' as text " posterior draws"
    display as text "CI level: `level'%"
end

* Drop existing definitions so we can redefine under matastrict
capture mata: mata drop _hs_extract_cates()
capture mata: mata drop _hs_cates_from_draws()

mata:
mata set matastrict on

// Note: _hs_quantile_pair() is provided by _hs_utils.ado (loaded above)


// -------------------------------------------------------------------------
// _hs_cates_from_draws(betavars, xtest_matname, xsd_matname,
//                      n_x, n_d, lo_prob, hi_prob)
//
// Bridge function called from the horseshoe_cates Stata command.
// Reads beta draws from current Stata data (after preserve/use),
// reads X_test and X_sd from named Stata matrices, and calls
// _hs_extract_cates(). Results are stored as Stata matrices
// (__hs_cate_hat, __hs_cate_lo, __hs_cate_hi).
// -------------------------------------------------------------------------
void _hs_cates_from_draws(
    string scalar betavars,
    string scalar xtest_matname,
    string scalar xsd_matname,
    real scalar   n_x,
    real scalar   n_d,
    real scalar   lo_prob,
    real scalar   hi_prob)
{
    real matrix    beta_all, X_test
    real rowvector X_sd

    // Read beta draws from current dataset (n_mcmc × p)
    beta_all = st_data(., tokens(betavars))

    // Read X_test and X_sd from Stata matrices
    X_test = st_matrix(xtest_matname)
    X_sd   = st_matrix(xsd_matname)

    // Call the core extraction function
    _hs_extract_cates(beta_all, X_test, X_sd, n_x, n_d, lo_prob, hi_prob)
}


// -------------------------------------------------------------------------
// _hs_extract_cates(beta_all, X_test, X_sd, n_x, n_d, lo_prob, hi_prob)
//
// Extract individual-level CATEs from posterior beta draws for a multi-arm
// treatment design. This is the shared implementation used by both the
// CATE wrapper (04_horseshoe_cate.do) and the validation script
// (04_horseshoe_test.do).
//
// Arguments:
//   beta_all  — n_mcmc × p matrix of posterior draws (SCALED parameterization)
//   X_test    — n_test × n_x matrix of test covariates (UNSCALED)
//   X_sd      — 1 × p rowvector of column SDs used for scaling
//   n_x       — number of covariate columns (54 for ACIC)
//   n_d       — number of treatment dummies (4 for ACIC)
//   lo_prob   — lower quantile probability (e.g., 0.025 for 95% CI)
//   hi_prob   — upper quantile probability (e.g., 0.975 for 95% CI)
//
// Returns (via Stata matrices, accessed with st_matrix()):
//     __hs_cate_hat — n_test × n_d: posterior mean iCATEs
//     __hs_cate_lo  — n_test × n_d: lower credible interval bounds
//     __hs_cate_hi  — n_test × n_d: upper credible interval bounds
//
//   Results are stored as Stata matrices (not Mata externals) so that
//   callers in mata { } blocks can read them with st_matrix("__hs_cate_hat").
//   This avoids a compile-time name resolution issue where Mata's { } block
//   compiler rejects references to external variables that don't yet exist.
//
// Column layout assumption (matching build_linear_data / 03_build_horseshoe_design):
//   Cols 1..n_x:          X covariate main effects
//   Cols n_x+1..n_x+n_d:  D treatment dummies (D_b, D_c, D_d, D_e)
//   Cols n_x+n_d+1..end:  X:D interactions, blocked by arm:
//     n_x cols for arm 1 (b), then n_x for arm 2 (c), etc.
// -------------------------------------------------------------------------
void _hs_extract_cates(
    real matrix    beta_all,
    real matrix    X_test,
    real rowvector X_sd,
    real scalar    n_x,
    real scalar    n_d,
    real scalar    lo_prob,
    real scalar    hi_prob)
{
    real scalar    n_mcmc, p, n_test, w_idx, i
    real scalar    gamma_col, delta_start, delta_end, j
    real matrix    beta_rescaled, delta_draws, tau_draws
    real matrix    cate_hat, cate_lo, cate_hi
    real colvector gamma_draws, draws_i
    real rowvector qi

    n_mcmc = rows(beta_all)
    p      = cols(beta_all)
    n_test = rows(X_test)

    // --- Rescale beta draws from scaled to original parameterization ---
    // beta_orig = beta_scaled / X_sd
    beta_rescaled = beta_all
    for (j = 1; j <= p; j++) {
        if (X_sd[j] != 0 & X_sd[j] != 1) {
            beta_rescaled[., j] = beta_rescaled[., j] / X_sd[j]
        }
    }

    // --- Allocate output matrices ---
    cate_hat = J(n_test, n_d, .)
    cate_lo  = J(n_test, n_d, .)
    cate_hi  = J(n_test, n_d, .)

    // --- Extract CATEs for each contrast ---
    for (w_idx = 1; w_idx <= n_d; w_idx++) {
        // gamma_w: treatment dummy coefficient
        gamma_col = n_x + w_idx

        // delta_w: interaction coefficients for this arm
        delta_start = n_x + n_d + (w_idx - 1) * n_x + 1
        delta_end   = delta_start + n_x - 1

        // gamma draws: n_mcmc × 1
        gamma_draws = beta_rescaled[., gamma_col]

        // delta draws: n_mcmc × n_x
        delta_draws = beta_rescaled[., (delta_start..delta_end)]

        // CATE draws: n_test × n_mcmc
        // tau_w^(s)(x_i) = gamma_w^(s) + x_i' * delta_w^(s)
        // X_test * delta_draws' gives n_test × n_mcmc
        // Then add gamma_w^(s) via row broadcast (:+ transposes 1×n_mcmc)
        tau_draws = X_test * delta_draws' :+ gamma_draws'

        // Posterior means: average across MCMC draws
        cate_hat[., w_idx] = mean(tau_draws')'

        // Credible intervals: quantiles per observation (single-sort)
        for (i = 1; i <= n_test; i++) {
            draws_i = tau_draws[i, .]'
            qi = _hs_quantile_pair(draws_i, lo_prob, hi_prob)
            cate_lo[i, w_idx] = qi[1]
            cate_hi[i, w_idx] = qi[2]
        }
    }

    // --- Store results as Stata matrices ---
    // Using st_matrix() instead of Mata external variables avoids the
    // compile-time name resolution issue in mata { } blocks. Callers
    // read these back with st_matrix("__hs_cate_hat"), etc.
    st_matrix("__hs_cate_hat", cate_hat)
    st_matrix("__hs_cate_lo",  cate_lo)
    st_matrix("__hs_cate_hi",  cate_hi)
}


end
