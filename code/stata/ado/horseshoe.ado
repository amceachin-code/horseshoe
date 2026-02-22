*! horseshoe.ado — Bayesian linear regression with the horseshoe prior
*! Version 1.1.0 — 2026-02-21
*!
*! Implements the auxiliary-variable Gibbs sampler from:
*!   Makalic & Schmidt (2015), "A simple sampler for the horseshoe estimator"
*!   arXiv:1508.03884v4
*!
*! Key feature: the penalized() option accepts a 1×p matrix of 0/1 flags.
*!   penalized[1,j] = 1  --> beta_j gets the horseshoe prior (shrinkage)
*!   penalized[1,j] = 0  --> beta_j gets a flat (improper) prior (no shrinkage)
*!
*! This lets you leave specific coefficients (e.g., treatment main effects)
*! unpenalized while applying the horseshoe prior to others.
*!
*! Author: Andrew McEachin
*! License: MIT


* ===========================================================================
* PROGRAM DEFINITION — horseshoe
*
* This is a Stata e-class estimation command. It:
*   1. Parses syntax (depvar, indepvars, options)
*   2. Converts data to Mata matrices
*   3. Calls the Mata Gibbs sampler (_hs_gibbs_run)
*   4. Posts results back to Stata's e() via _hs_post_results
*   5. Optionally saves full posterior draws to a .dta file
* ===========================================================================

program define horseshoe, eclass
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

    * -------------------------------------------------------------------
    * REPLAY: if called with no arguments, re-display last results
    * -------------------------------------------------------------------
    if replay() {
        if "`e(cmd)'" != "horseshoe" {
            error 301
        }
        _hs_display `0'
        exit
    }

    * -------------------------------------------------------------------
    * SYNTAX PARSING
    *
    * Syntax:
    *   horseshoe depvar indepvars [if] [in], [options]
    *
    * The `varlist(min=2)` requires at least 2 variables: 1 depvar + 1+ indepvars.
    * -------------------------------------------------------------------
    syntax varlist(min=2 numeric) [if] [in], ///
        [                                     ///
        PENalized(name)                       /// 1×p Stata matrix of 0/1 penalization flags
        NMCMC(integer 1000)                   /// Number of post-burnin MCMC draws to store
        BURNin(integer 500)                   /// Number of burnin iterations to discard
        THIN(integer 1)                       /// Thinning interval (keep every thin-th draw)
        LAMBDAScale(real 1)                   /// Half-Cauchy scale on local shrinkage lambda_j
        TAUScale(real 1)                      /// Half-Cauchy scale on global shrinkage tau
        SAVing(string)                        /// File path to save full posterior draws (.dta)
        REPLACE                               /// Allow overwriting the saving() file
        SEED(string)                          /// RNG seed for reproducibility
        Level(cilevel)                        /// Confidence/credible interval level (default 95)
        VERBose                               /// Print iteration progress every 500 iterations
        ]

    * -------------------------------------------------------------------
    * INPUT VALIDATION
    * -------------------------------------------------------------------

    * Separate depvar from indepvars
    gettoken depvar indepvars : varlist
    local p : word count `indepvars'

    if `p' < 1 {
        display as error "At least one independent variable is required."
        exit 198
    }

    * Validate MCMC parameters
    if `nmcmc' < 1 {
        display as error "nmcmc() must be a positive integer, got `nmcmc'"
        exit 198
    }
    if `burnin' < 0 {
        display as error "burnin() must be non-negative, got `burnin'"
        exit 198
    }
    if `thin' < 1 {
        display as error "thin() must be a positive integer, got `thin'"
        exit 198
    }

    * Validate scale parameters
    if `lambdascale' <= 0 {
        display as error "lambdascale() must be positive, got `lambdascale'"
        exit 198
    }
    if `tauscale' <= 0 {
        display as error "tauscale() must be positive, got `tauscale'"
        exit 198
    }

    * Validate penalized matrix dimensions (if provided)
    if "`penalized'" != "" {
        * Check the matrix exists
        capture confirm matrix `penalized'
        if _rc {
            display as error "Matrix `penalized' not found."
            exit 111
        }
        * Check it's 1×p
        if rowsof(`penalized') != 1 | colsof(`penalized') != `p' {
            display as error ///
                "penalized() matrix must be 1 × `p', " ///
                "got " rowsof(`penalized') " × " colsof(`penalized')
            exit 503
        }
        * Check that all entries are exactly 0 or 1
        tempname pen_check
        matrix `pen_check' = `penalized'
        forvalues j = 1/`p' {
            if `pen_check'[1, `j'] != 0 & `pen_check'[1, `j'] != 1 {
                display as error ///
                    "penalized() matrix must contain only 0s and 1s; " ///
                    "found " `pen_check'[1, `j'] " at position `j'"
                exit 198
            }
        }
    }

    * Validate saving() option
    if "`saving'" != "" & "`replace'" == "" {
        * Check if file already exists
        capture confirm file "`saving'"
        if !_rc {
            display as error ///
                "File `saving' already exists. " ///
                "Specify replace to overwrite."
            exit 602
        }
    }

    * -------------------------------------------------------------------
    * SET RNG SEED (if requested)
    * -------------------------------------------------------------------
    if "`seed'" != "" {
        * Validate seed: must be a non-negative integer or a valid state string
        capture confirm integer number `seed'
        if _rc == 0 {
            if `seed' < 0 {
                display as error "seed() must be a non-negative integer, got `seed'"
                exit 198
            }
            set seed `seed'
        }
        else {
            * Try as a rng_state string (from e(rng_state))
            capture set rngstate `seed'
            if _rc {
                display as error ///
                    "seed() must be a non-negative integer or a valid RNG state string"
                exit 198
            }
        }
    }

    * Capture current RNG state for storage in e()
    local rng_current_state = c(rngstate)

    * -------------------------------------------------------------------
    * MARK ESTIMATION SAMPLE
    * -------------------------------------------------------------------
    marksample touse
    quietly count if `touse'
    local n = r(N)

    if `n' < 2 {
        display as error "Insufficient observations (n = `n')."
        exit 2001
    }

    * Check propriety: count unpenalized parameters
    if "`penalized'" != "" {
        * Count zeros in penalized matrix
        tempname pen_mat
        matrix `pen_mat' = `penalized'
        local p_free = 0
        forvalues j = 1/`p' {
            if `pen_mat'[1, `j'] == 0 {
                local p_free = `p_free' + 1
            }
        }
        if `p_free' > 0 & `n' <= `p_free' {
            display as text ///
                "Warning: n (`n') <= number of unpenalized parameters (`p_free')." _newline ///
                "  Posterior may be improper for unpenalized coefficients."
        }
    }

    * -------------------------------------------------------------------
    * TRANSFER DATA TO MATA
    *
    * We pass the depvar, indepvars, touse marker, penalized matrix,
    * and all scalar options to the Mata Gibbs sampler.
    * -------------------------------------------------------------------

    * Store variable names for later labeling
    local varnames `indepvars'

    * Compute CI bounds from the level option
    * level is e.g. 95 → alpha = 0.05 → lower = 0.025, upper = 0.975
    local alpha = (100 - `level') / 100
    local lower_pct = `alpha' / 2
    local upper_pct = 1 - `alpha' / 2

    * -------------------------------------------------------------------
    * CALL MATA GIBBS SAMPLER
    *
    * The Mata function _hs_gibbs_run() does all the heavy lifting.
    * It reads data from Stata via st_data(), runs the Gibbs sampler,
    * and returns a struct with all posterior draws.
    *
    * Then _hs_post_results() pushes the summarized results back to
    * Stata's e() class.
    * -------------------------------------------------------------------

    * Set verbose flag for Mata
    if "`verbose'" != "" {
        local verbose_flag = 1
    }
    else {
        local verbose_flag = 0
    }

    * Set penalized matrix name for Mata (empty string if not specified)
    local pen_matname `penalized'

    * Call the main Mata entry point
    * Wrap in capture noisily so we can clean up __hs_* globals on error
    capture noisily mata: _hs_stata_main( ///
        "`depvar'",                      ///
        "`indepvars'",                   ///
        "`touse'",                       ///
        "`pen_matname'",                 ///
        `nmcmc',                         ///
        `burnin',                        ///
        `thin',                          ///
        `lambdascale',                   ///
        `tauscale',                      ///
        `lower_pct',                     ///
        `upper_pct',                     ///
        `verbose_flag',                  ///
        "`saving'",                      ///
        "`varnames'"                     ///
    )
    if _rc {
        * Clean up any stale __hs_* globals before re-raising the error
        capture matrix drop __hs_b __hs_b_sd __hs_b_lower __hs_b_upper __hs_V
        capture scalar drop __hs_sigma2 __hs_tau2 __hs_N __hs_p __hs_p_pen
        capture scalar drop __hs_nmcmc __hs_lambda_scale __hs_tau_scale
        capture scalar drop __hs_burnin __hs_thin __hs_n_chol_fallback
        capture scalar drop __hs_ess_min __hs_ess_sigma2 __hs_ess_tau2
        capture matrix drop __hs_ess_beta
        exit _rc
    }

    * -------------------------------------------------------------------
    * POST RESULTS TO E-CLASS
    *
    * Mata stored results in regular matrices/scalars (__hs_*).
    * Now move them into the e() class via ereturn commands.
    * -------------------------------------------------------------------

    * Initialize e() with the coefficient vector and variance-covariance matrix
    * ereturn post requires e(b); posting e(V) simultaneously unlocks
    * test, lincom, nlcom, margins, and predict stdp.
    tempname b_post V_post
    matrix `b_post' = __hs_b
    matrix `V_post' = __hs_V

    * Label the coefficient and variance matrices with variable names
    matrix colnames `b_post' = `indepvars'
    matrix colnames `V_post' = `indepvars'
    matrix rownames `V_post' = `indepvars'

    * Post the coefficient vector and variance matrix to initialize e()
    ereturn post `b_post' `V_post', esample(`touse')

    * Now post additional matrices
    tempname sd_post lo_post hi_post
    matrix `sd_post' = __hs_b_sd
    matrix `lo_post' = __hs_b_lower
    matrix `hi_post' = __hs_b_upper
    ereturn matrix b_sd    = `sd_post'
    ereturn matrix b_lower = `lo_post'
    ereturn matrix b_upper = `hi_post'

    * Store penalization flags for replay display
    tempname pen_post
    if "`penalized'" != "" {
        matrix `pen_post' = `penalized'
    }
    else {
        matrix `pen_post' = J(1, `p', 1)
    }
    matrix colnames `pen_post' = `indepvars'
    ereturn matrix penalized = `pen_post'

    * Post scalars
    ereturn scalar sigma2       = scalar(__hs_sigma2)
    ereturn scalar tau2         = scalar(__hs_tau2)
    ereturn scalar N            = scalar(__hs_N)
    ereturn scalar p            = scalar(__hs_p)
    ereturn scalar p_pen        = scalar(__hs_p_pen)
    ereturn scalar nmcmc        = scalar(__hs_nmcmc)
    ereturn scalar lambda_scale = scalar(__hs_lambda_scale)
    ereturn scalar tau_scale    = scalar(__hs_tau_scale)
    ereturn scalar burnin         = scalar(__hs_burnin)
    ereturn scalar thin           = scalar(__hs_thin)
    ereturn scalar n_chol_fallback = scalar(__hs_n_chol_fallback)

    * Post convergence diagnostics
    ereturn scalar ess_min    = scalar(__hs_ess_min)
    ereturn scalar ess_sigma2 = scalar(__hs_ess_sigma2)
    ereturn scalar ess_tau2   = scalar(__hs_ess_tau2)
    tempname ess_beta_post
    matrix `ess_beta_post' = __hs_ess_beta
    matrix colnames `ess_beta_post' = `indepvars'
    ereturn matrix ess_beta  = `ess_beta_post'

    * Post strings
    ereturn local cmd       "horseshoe"
    ereturn local cmdline   `"horseshoe `0'"'
    ereturn local predict   "horseshoe_p"
    ereturn local depvar    "`__hs_depvar'"
    ereturn local seed      "`seed'"
    ereturn local rng_state "`rng_current_state'"

    * Clean up temporary Stata matrices and scalars
    capture matrix drop __hs_b __hs_b_sd __hs_b_lower __hs_b_upper __hs_V
    capture scalar drop __hs_sigma2 __hs_tau2 __hs_N __hs_p __hs_p_pen
    capture scalar drop __hs_nmcmc __hs_lambda_scale __hs_tau_scale
    capture scalar drop __hs_burnin __hs_thin __hs_n_chol_fallback
    capture scalar drop __hs_ess_min __hs_ess_sigma2 __hs_ess_tau2
    capture matrix drop __hs_ess_beta

    * -------------------------------------------------------------------
    * DISPLAY RESULTS TABLE
    * -------------------------------------------------------------------
    _hs_display, level(`level')

    * Report saving file
    if "`saving'" != "" {
        display ""
        display as text "Posterior draws saved to: " as result "`saving'"
    }

end


* ===========================================================================
* _hs_display — Subroutine for displaying horseshoe results
*
* Reads everything from e() so it works for both initial display and replay.
* ===========================================================================
program define _hs_display
    syntax [, Level(cilevel)]

    display ""
    display as text "Horseshoe prior regression" _col(56) ///
        "Number of obs   = " as result %8.0f e(N)
    display as text _col(56) ///
        "Predictors      = " as result %8.0f e(p)
    display as text _col(56) ///
        "Penalized       = " as result %8.0f e(p_pen)
    display as text _col(56) ///
        "MCMC draws      = " as result %8.0f e(nmcmc)
    display ""
    display as text "Posterior mean of sigma^2 = " ///
        as result %9.4f e(sigma2)
    display as text "Posterior mean of tau^2   = " ///
        as result %9.6f e(tau2)
    display ""
    display as text "Min effective sample size = " ///
        as result %9.0f e(ess_min)
    display as text "ESS for sigma^2           = " ///
        as result %9.0f e(ess_sigma2)
    display as text "ESS for tau^2             = " ///
        as result %9.0f e(ess_tau2)

    * Warn if ESS is suspiciously low
    if e(ess_min) < 100 {
        display as error ///
            "Warning: minimum ESS < 100. " ///
            "Consider increasing nmcmc() or burnin()."
    }
    display ""

    * Display coefficient table
    * Header
    display as text "{hline 13}{c +}{hline 64}"
    display as text %12s "Variable" " {c |}" ///
        %12s "Post.Mean" ///
        %12s "Post.SD" ///
        %12s "Lower CI" ///
        %12s "Upper CI" ///
        %12s "Penalized"
    display as text "{hline 13}{c +}{hline 64}"

    * Body — loop over variables using column names from e(b)
    tempname b_vec sd_vec lo_vec hi_vec pen_row
    matrix `b_vec'   = e(b)
    matrix `sd_vec'  = e(b_sd)
    matrix `lo_vec'  = e(b_lower)
    matrix `hi_vec'  = e(b_upper)
    matrix `pen_row' = e(penalized)

    local varnames : colnames `b_vec'

    local i = 0
    foreach var of local varnames {
        local ++i
        local pen_str "Yes"
        if `pen_row'[1, `i'] == 0 {
            local pen_str "No"
        }
        display as text %12s abbrev("`var'", 12) " {c |}" ///
            as result ///
            %12.4f `b_vec'[1, `i'] ///
            %12.4f `sd_vec'[1, `i'] ///
            %12.4f `lo_vec'[1, `i'] ///
            %12.4f `hi_vec'[1, `i'] ///
            as text %12s "`pen_str'"
    }

    display as text "{hline 13}{c +}{hline 64}"
    display as text "CI level: `level'%"
end


* ===========================================================================
* MATA CODE — Embedded functions for the Gibbs sampler
*
* Functions:
*   _hs_stata_main()     — Entry point called from Stata; orchestrates everything
*   _hs_gibbs_run()      — Core Gibbs sampler loop (6 conditional steps)
*   _hs_inv_gamma()      — Single draw from InverseGamma(shape, rate)
*   _hs_ess()            — Effective sample size (initial positive sequence)
*   _hs_convergence_diag() — Compute ESS for all parameters
*   _hs_post_results()   — Push Mata struct results → Stata matrices
*   _hs_save_draws()     — Save posterior draws to .dta file
*
* External (from _hs_utils.ado):
*   _hs_quantile_pair()  — Paired lo/hi quantile in a single sort pass
*
* Struct:
*   _hs_result           — Holds all posterior draws and metadata
* ===========================================================================

mata:
mata set matastrict on


// -------------------------------------------------------------------------
// Struct definition: holds all posterior draws and sampler metadata
// -------------------------------------------------------------------------
struct _hs_result {
    // Posterior draws
    real matrix    beta_draws      // p × n_mcmc: posterior draws of beta
    real colvector sigma2_draws    // n_mcmc × 1: posterior draws of sigma^2
    real colvector tau2_draws      // n_mcmc × 1: posterior draws of tau^2
    real matrix    lambda2_draws   // p_pen × n_mcmc: posterior draws of lambda_j^2

    // Index mapping: which columns of X are penalized
    real colvector pen_idx         // p_pen × 1: indices of penalized columns (1-based)
    real colvector free_idx        // p_free × 1: indices of unpenalized columns

    // Dimensions
    real scalar    n               // number of observations
    real scalar    p               // total number of predictors
    real scalar    p_pen           // number of penalized predictors
    real scalar    n_mcmc          // number of stored MCMC draws
}


// -------------------------------------------------------------------------
// _hs_inv_gamma(shape, rate)
//
// Draw a single scalar from InverseGamma(shape, rate).
//
// Parameterization:
//   If X ~ InvGamma(a, b), then 1/X ~ Gamma(a, 1/b) where b is the RATE.
//   Mata's rgamma(rows, cols, shape, scale) returns an r×c matrix of draws
//   from Gamma(shape, scale).
//   So: 1/X ~ Gamma(shape, 1/rate)  →  rgamma(1, 1, shape, 1/rate)
//   Then X = 1/rgamma(1, 1, shape, 1/rate).
//
// This matches R's: 1/rgamma(1, shape=shape, rate=rate)
// Since R's rgamma uses shape+rate, and rgamma(1, shape, rate=r)
// = rgamma(1, shape, scale=1/r).
// -------------------------------------------------------------------------
real scalar _hs_inv_gamma(real scalar shape, real scalar rate)
{
    real scalar x

    // rgamma(rows, cols, shape, scale): draw 1×1 scalar
    x = rgamma(1, 1, shape, 1 / rate)

    // Guard against zero draws (numerically possible with very large rates)
    if (x < 1e-300) x = 1e-300

    return(1 / x)
}



// -------------------------------------------------------------------------
// _hs_gibbs_run(y, X, pen_flag, n_mcmc, burnin, thin,
//               lambda_scale, tau_scale, verbose)
//
// Core Gibbs sampler implementing Makalic & Schmidt (2015), Algorithm 1,
// with the extension that only a subset of coefficients receive the
// horseshoe prior (pen_flag[j] = 1). Unpenalized coefficients get a
// flat (improper) prior.
//
// Arguments:
//   y            — n × 1 outcome vector (should be centered if no intercept)
//   X            — n × p design matrix
//   pen_flag     — 1 × p row vector: 1 = penalized, 0 = unpenalized
//   n_mcmc       — number of post-burnin draws to store
//   burnin       — number of initial iterations to discard
//   thin         — thinning interval
//   lambda_scale — scale of half-Cauchy prior on lambda_j
//   tau_scale    — scale of half-Cauchy prior on tau
//   verbose      — 1 = print progress every 500 iterations
//
// Returns:
//   struct _hs_result with all posterior draws
// -------------------------------------------------------------------------
struct _hs_result scalar _hs_gibbs_run(
    real colvector y,
    real matrix    X,
    real rowvector pen_flag,
    real scalar    n_mcmc,
    real scalar    burnin,
    real scalar    thin,
    real scalar    lambda_scale,
    real scalar    tau_scale,
    real scalar    verbose)
{
    struct _hs_result scalar res

    // --- Dimensions ---
    real scalar n, p, p_pen
    n = rows(X)
    p = cols(X)

    // --- Build penalized / free index vectors ---
    // pen_idx: column indices where pen_flag == 1
    // free_idx: column indices where pen_flag == 0
    real colvector pen_idx, free_idx
    real scalar j, n_pen, n_free

    // Count penalized and free
    n_pen = 0
    n_free = 0
    for (j = 1; j <= p; j++) {
        if (pen_flag[j] == 1) n_pen++
        else n_free++
    }
    p_pen = n_pen

    if (p_pen == 0) {
        errprintf("No penalized parameters. Use OLS instead.\n")
        exit(198)
    }

    pen_idx  = J(p_pen, 1, .)
    if (n_free > 0) {
        free_idx = J(n_free, 1, .)
    }
    else {
        free_idx = J(0, 1, .)
    }

    n_pen = 0
    n_free = 0
    for (j = 1; j <= p; j++) {
        if (pen_flag[j] == 1) {
            pen_idx[++n_pen] = j
        }
        else {
            free_idx[++n_free] = j
        }
    }

    // --- Precompute fixed quantities ---
    real matrix    XtX
    real colvector Xty

    XtX = cross(X, X)     // p × p: X'X
    Xty = cross(X, y)     // p × 1: X'y

    // Half-Cauchy scale-mixture constants for auxiliary variable conditionals:
    //   nu_j | . ~ IG(1, 1/lambda_j^2 + 1/lambda_scale^2)
    //   xi   | . ~ IG(1, 1/tau^2       + 1/tau_scale^2)
    real scalar inv_lambda_scale_sq, inv_tau_scale_sq
    inv_lambda_scale_sq = 1 / (lambda_scale^2)
    inv_tau_scale_sq    = 1 / (tau_scale^2)

    // Total iterations: burnin + (n_mcmc * thin)
    real scalar total_iter
    total_iter = burnin + n_mcmc * thin

    // --- Initialize parameters ---

    // Beta: ridge estimate (lambda=1) for numerical stability at iteration 1
    // Solves (X'X + I)beta = X'y
    real colvector beta
    real matrix    A_init
    A_init = XtX + I(p)
    beta = cholsolve(A_init, Xty)
    // Fallback: if cholsolve fails (singular), use lusolve
    if (beta[1] == .) {
        beta = lusolve(A_init, Xty)
    }

    // Sigma^2: initialize from residuals
    real scalar sigma2
    real colvector resid_init
    resid_init = y - X * beta
    sigma2 = cross(resid_init, resid_init) / max((n - p, 1))

    // Local shrinkage: lambda_j^2 for penalized params (all start at 1)
    real colvector lambda2
    lambda2 = J(p_pen, 1, 1)

    // Global shrinkage: tau^2
    real scalar tau2
    tau2 = 1

    // Auxiliary variables for half-Cauchy decomposition
    real colvector nu       // one per penalized param
    real scalar    xi
    nu = J(p_pen, 1, 1)
    xi = 1

    // --- Allocate storage for posterior draws ---
    real matrix    beta_store
    real colvector sigma2_store, tau2_store
    real matrix    lambda2_store

    beta_store    = J(p, n_mcmc, .)
    sigma2_store  = J(n_mcmc, 1, .)
    tau2_store    = J(n_mcmc, 1, .)
    lambda2_store = J(p_pen, n_mcmc, .)

    real scalar save_idx
    save_idx = 0

    // --- Scratch variables for the Gibbs loop ---
    real colvector diag_penalty   // p × 1: diagonal penalty vector
    real matrix    A               // p × p: precision matrix X'X + diag(penalty)
    real matrix    L               // p × p: lower Cholesky factor of A
    real colvector mu_beta         // p × 1: posterior mean of beta
    real colvector z_draw          // p × 1: standard normal draw
    real colvector resid           // n × 1: residual y - Xb
    real scalar    rss             // residual sum of squares
    real scalar    pen_term        // sum of penalized beta^2 / (tau2*lambda2)
    real scalar    shape_sigma, rate_sigma
    real scalar    shape_tau, rate_tau
    real colvector beta_pen           // p_pen × 1: penalized beta subvector
    real colvector lambda_rates       // p_pen × 1: vectorized rates for lambda^2
    real colvector nu_rates           // p_pen × 1: vectorized rates for nu
    real scalar    iter, jj
    real scalar    n_chol_fallback   // count Cholesky fallback events

    // =======================================================================
    // GIBBS SAMPLER MAIN LOOP
    //
    // 6 conditional steps per iteration, matching the R implementation exactly.
    // See horseshoe_gibbs.R for detailed derivations of each conditional.
    // =======================================================================

    n_chol_fallback = 0

    if (verbose) {
        printf("{txt}  [Horseshoe Gibbs] Starting sampler: %g total iterations ", total_iter)
        printf("(%g burnin + %g × %g thin)\n", burnin, n_mcmc, thin)
        printf("{txt}  [Horseshoe Gibbs] n = %g, p = %g (%g penalized, %g free)\n",
               n, p, p_pen, rows(free_idx))
        printf("{txt}  [Horseshoe Gibbs] lambda_scale = %9.2f, tau_scale = %9.4f\n",
               lambda_scale, tau_scale)
        displayflush()
    }

    for (iter = 1; iter <= total_iter; iter++) {

        // -----------------------------------------------------------------
        // Step 1: Sample beta | rest
        //
        // beta | . ~ N_p(A^{-1}X'y, sigma^2 * A^{-1})
        // where A = X'X + diag(penalty)
        //
        // Penalty diagonal:
        //   1/(tau^2 * lambda_j^2) for penalized j
        //   0                      for unpenalized j (flat prior)
        //
        // Sampling via Cholesky:
        //   L = cholesky(A)     [lower triangular: A = LL']
        //   mu = L'\(L\X'y)    [reuse L; avoids re-factoring A]
        //   z ~ N(0, I_p)
        //   beta = mu + sqrt(sigma2) * solveupper(L', z)
        //
        // Note: Mata cholesky() returns LOWER triangular (unlike R's chol()
        //   which returns upper). lusolve(L', z) solves L'x = z, i.e.,
        //   x = (L')^{-1} z = (L^{-1})' z. Since A = LL',
        //   A^{-1} = (L')^{-1} L^{-1}, so
        //   sqrt(sigma2) * (L')^{-1} z has covariance sigma2 * A^{-1}. ✓
        // -----------------------------------------------------------------
        diag_penalty = J(p, 1, 0)
        diag_penalty[pen_idx] = 1 :/ (tau2 * lambda2)

        A = XtX + diag(diag_penalty)
        L = cholesky(A)

        // If Cholesky fails (shouldn't with ridge penalty, but be safe),
        // add a small ridge and retry
        if (L[1,1] == .) {
            A = A + 1e-8 * I(p)
            L = cholesky(A)
            n_chol_fallback++
            if (L[1,1] == .) {
                errprintf("Cholesky factorization failed even with ridge adjustment at iteration %g.\n", iter)
                errprintf("  The design matrix may be severely rank-deficient.\n")
                exit(error(430))
            }
        }

        mu_beta = solveupper(L', solvelower(L, Xty))
        z_draw = rnormal(p, 1, 0, 1)
        beta = mu_beta + sqrt(sigma2) * solveupper(L', z_draw)

        // -----------------------------------------------------------------
        // Step 2: Sample sigma^2 | rest
        //
        // sigma^2 | . ~ IG((n + p_pen)/2,
        //                   [||y - Xb||^2 + sum_pen beta_j^2/(tau^2*lambda_j^2)] / 2)
        //
        // Only penalized betas contribute the prior penalty term.
        // -----------------------------------------------------------------
        resid = y - X * beta
        rss = cross(resid, resid)

        beta_pen = beta[pen_idx]
        pen_term = sum(beta_pen :^2 :/ (tau2 * lambda2))

        shape_sigma = (n + p_pen) / 2
        rate_sigma  = (rss + pen_term) / 2
        sigma2 = _hs_inv_gamma(shape_sigma, rate_sigma)

        // Numerical safety: check sigma2 for NaN/Inf
        if (missing(sigma2) | sigma2 > 1e200) {
            errprintf("Gibbs sampler diverged at iteration %g: sigma2 = %g\n", iter, sigma2)
            errprintf("  This usually indicates near-collinear design columns or\n")
            errprintf("  a highly ill-conditioned X'X matrix.\n")
            exit(error(430))
        }

        // -----------------------------------------------------------------
        // Step 3: Sample lambda_j^2 | rest  (loop over penalized j)
        //
        // lambda_j^2 | . ~ IG(1, 1/nu_j + beta_j^2 / (2*tau^2*sigma^2))
        //
        // Must loop because each lambda_j has a different rate parameter.
        // Mata rgamma() is scalar, so vectorization isn't possible here.
        // -----------------------------------------------------------------
        lambda_rates = 1 :/ nu + beta_pen :^2 :/ (2 * tau2 * sigma2)
        for (jj = 1; jj <= p_pen; jj++) {
            lambda2[jj] = _hs_inv_gamma(1, lambda_rates[jj])
        }

        // -----------------------------------------------------------------
        // Step 4: Sample tau^2 | rest
        //
        // tau^2 | . ~ IG((p_pen+1)/2,
        //                 1/xi + sum_pen(beta_j^2 / lambda_j^2) / (2*sigma^2))
        //
        // The sum runs only over penalized coefficients.
        // -----------------------------------------------------------------
        shape_tau = (p_pen + 1) / 2
        // rate = 1/xi + sum_pen(beta_j^2 / lambda_j^2) / (2*sigma^2)
        // Note: the 1/xi term is NOT divided by (2*sigma^2)
        pen_term = sum(beta_pen :^2 :/ lambda2)
        rate_tau = 1 / xi + pen_term / (2 * sigma2)
        tau2 = _hs_inv_gamma(shape_tau, rate_tau)

        // Numerical safety: check tau2 for NaN/Inf
        if (missing(tau2) | tau2 > 1e200) {
            errprintf("Gibbs sampler diverged at iteration %g: tau2 = %g\n", iter, tau2)
            errprintf("  This usually indicates a degenerate posterior for the\n")
            errprintf("  global shrinkage parameter.\n")
            exit(error(430))
        }

        // -----------------------------------------------------------------
        // Step 5: Sample nu_j | rest  (auxiliary for lambda_j)
        //
        // nu_j | . ~ IG(1, 1/lambda_scale^2 + 1/lambda_j^2)
        //
        // These are the auxiliary variables from the half-Cauchy
        // scale-mixture representation:
        //   lambda_j ~ C+(0, lambda_scale) iff
        //   lambda_j^2 | nu_j ~ IG(1/2, 1/nu_j) and
        //   nu_j ~ IG(1/2, 1/lambda_scale^2)
        //
        // The combined conditional is IG(1, 1/lambda_scale^2 + 1/lambda_j^2).
        // -----------------------------------------------------------------
        nu_rates = J(p_pen, 1, inv_lambda_scale_sq) + 1 :/ lambda2
        for (jj = 1; jj <= p_pen; jj++) {
            nu[jj] = _hs_inv_gamma(1, nu_rates[jj])
        }

        // -----------------------------------------------------------------
        // Step 6: Sample xi | rest  (auxiliary for tau)
        //
        // xi | . ~ IG(1, 1/tau_scale^2 + 1/tau^2)
        //
        // Same half-Cauchy decomposition as nu_j, but for the global
        // shrinkage parameter tau.
        // -----------------------------------------------------------------
        xi = _hs_inv_gamma(1, inv_tau_scale_sq + 1 / tau2)

        // -----------------------------------------------------------------
        // Store draws (after burnin, respecting thinning)
        //
        // A draw is stored when:
        //   iter > burnin AND (iter - burnin) mod thin == 0
        // -----------------------------------------------------------------
        if (iter > burnin & mod(iter - burnin, thin) == 0) {
            save_idx++
            beta_store[., save_idx]    = beta
            sigma2_store[save_idx]     = sigma2
            tau2_store[save_idx]       = tau2
            lambda2_store[., save_idx] = lambda2
        }

        // -----------------------------------------------------------------
        // Progress reporting every 500 iterations
        // -----------------------------------------------------------------
        if (verbose & mod(iter, 500) == 0) {
            printf("{txt}  [Horseshoe Gibbs] iter %g / %g | tau = %9.4f | sigma = %9.4f\n",
                   iter, total_iter, sqrt(tau2), sqrt(sigma2))
            displayflush()
        }
    }

    // --- Post-loop warnings ---
    if (n_chol_fallback > 0 & verbose) {
        printf("{err}  Warning: Cholesky factorization required ridge adjustment in %g of %g iterations.\n",
               n_chol_fallback, total_iter)
        printf("{err}  This may indicate near-collinear columns in the design matrix.\n")
        displayflush()
    }

    // Check for extreme lambda2 values (indicates numerical issues)
    if (max(lambda2) > 1e100) {
        if (verbose) {
            printf("{err}  Warning: some lambda^2 values are extremely large (max = %g).\n", max(lambda2))
            printf("{err}  This may indicate problematic scaling of design columns.\n")
            displayflush()
        }
    }

    // Store Cholesky fallback count for e() posting
    st_numscalar("__hs_n_chol_fallback", n_chol_fallback)

    // --- Pack results into struct ---
    res.beta_draws    = beta_store
    res.sigma2_draws  = sigma2_store
    res.tau2_draws    = tau2_store
    res.lambda2_draws = lambda2_store
    res.pen_idx       = pen_idx
    res.free_idx      = free_idx
    res.n             = n
    res.p             = p
    res.p_pen         = p_pen
    res.n_mcmc        = n_mcmc

    return(res)
}


// -------------------------------------------------------------------------
// _hs_post_results(res, lower_pct, upper_pct, lambda_scale, tau_scale,
//                  burnin, thin, varnames)
//
// Compute posterior summaries from the Gibbs sampler output and post them
// to Stata's e() class as matrices and scalars.
//
// Stored results:
//   e(b)           — 1 × p: posterior means of beta
//   e(b_sd)        — 1 × p: posterior SDs of beta
//   e(b_lower)     — 1 × p: lower credible interval bound
//   e(b_upper)     — 1 × p: upper credible interval bound
//   e(sigma2)      — scalar: posterior mean of sigma^2
//   e(tau2)        — scalar: posterior mean of tau^2
//   e(N)           — scalar: number of observations
//   e(p)           — scalar: number of predictors
//   e(p_pen)       — scalar: number of penalized predictors
//   e(nmcmc)       — scalar: number of stored MCMC draws
//   e(lambda_scale)— scalar: lambda_scale used
//   e(tau_scale)   — scalar: tau_scale used
// -------------------------------------------------------------------------
void _hs_post_results(
    struct _hs_result scalar res,
    real scalar lower_pct,
    real scalar upper_pct,
    real scalar lambda_scale,
    real scalar tau_scale,
    real scalar burnin,
    real scalar thin,
    string scalar depvar,
    string scalar varnames)
{
    real scalar    p, n_mcmc, j
    real rowvector b_mean, b_sd, b_lower, b_upper, qi
    real colvector beta_j_draws

    p      = res.p
    n_mcmc = res.n_mcmc

    // --- Compute posterior summaries for each beta_j ---
    b_mean  = J(1, p, 0)
    b_sd    = J(1, p, 0)
    b_lower = J(1, p, 0)
    b_upper = J(1, p, 0)

    for (j = 1; j <= p; j++) {
        beta_j_draws = res.beta_draws[j, .]'  // n_mcmc × 1

        b_mean[j]  = mean(beta_j_draws)
        b_sd[j]    = sqrt(variance(beta_j_draws))
        qi = _hs_quantile_pair(beta_j_draws, lower_pct, upper_pct)
        b_lower[j] = qi[1]
        b_upper[j] = qi[2]
    }

    // --- Store matrices as regular Stata matrices (not e-class) ---
    // The Stata program will then use ereturn post/matrix to move them
    // into e(). Direct st_matrix("e(b)", ...) requires e() to be
    // initialized first, which isn't the case from Mata.
    st_matrix("__hs_b",       b_mean)
    st_matrix("__hs_b_sd",    b_sd)
    st_matrix("__hs_b_lower", b_lower)
    st_matrix("__hs_b_upper", b_upper)

    // --- Variance-covariance matrix of posterior draws ---
    // V = Var(beta_draws') is p × p, used for e(V) so that Stata's
    // test, lincom, nlcom, margins, and predict stdp all work.
    real matrix V
    V = variance(res.beta_draws')
    st_matrix("__hs_V", V)

    // --- Store scalars ---
    st_numscalar("__hs_sigma2",       mean(res.sigma2_draws))
    st_numscalar("__hs_tau2",         mean(res.tau2_draws))
    st_numscalar("__hs_N",            res.n)
    st_numscalar("__hs_p",            res.p)
    st_numscalar("__hs_p_pen",        res.p_pen)
    st_numscalar("__hs_nmcmc",        res.n_mcmc)
    st_numscalar("__hs_lambda_scale", lambda_scale)
    st_numscalar("__hs_tau_scale",    tau_scale)
    st_numscalar("__hs_burnin",       burnin)
    st_numscalar("__hs_thin",         thin)

    // Store strings in Stata locals for later ereturn global
    st_local("__hs_depvar", depvar)
    st_local("__hs_varnames", varnames)
}


// -------------------------------------------------------------------------
// _hs_save_draws(res, filename, varnames)
//
// Save the full posterior draws to a .dta file. The file contains:
//   draw     — draw number (1 to n_mcmc)
//   sigma2   — posterior draw of sigma^2
//   tau2     — posterior draw of tau^2
//   beta_1 ... beta_p — posterior draws of each coefficient
//
// This allows users to do custom posterior analysis (e.g., functions of
// parameters, predictive distributions, etc.).
// -------------------------------------------------------------------------
void _hs_save_draws(
    struct _hs_result scalar res,
    string scalar filename,
    string scalar varnames)
{
    real scalar n_mcmc, p, j
    string rowvector vnames

    n_mcmc = res.n_mcmc
    p      = res.p
    vnames = tokens(varnames)

    // Build the data matrix: n_mcmc × (3 + p)
    // Columns: draw, sigma2, tau2, beta_1, ..., beta_p
    real matrix data
    data = J(n_mcmc, 3 + p, .)

    // Column 1: draw number (1 to n_mcmc)
    data[., 1] = (1::n_mcmc)

    // Column 2: sigma2 draws
    data[., 2] = res.sigma2_draws

    // Column 3: tau2 draws
    data[., 3] = res.tau2_draws

    // Columns 4 to 3+p: beta draws (transposed from p × n_mcmc to n_mcmc × p)
    data[., (4..3+p)] = res.beta_draws'

    // --- Build variable name list ---
    // Use actual variable names (prefixed with "b_") for readability.
    // Sanitize: replace non-alphanumeric chars with underscore, truncate to
    // 30 chars (Stata 32-char limit minus "b_" prefix).
    string colvector all_names
    string scalar    nm
    all_names = J(3 + p, 1, "")
    all_names[1] = "draw"
    all_names[2] = "sigma2"
    all_names[3] = "tau2"
    for (j = 1; j <= p; j++) {
        if (j <= length(vnames)) {
            nm = vnames[j]
            // Sanitize: replace non-alphanumeric with underscore
            nm = subinstr(nm, ".", "_")
            nm = subinstr(nm, ":", "_")
            nm = subinstr(nm, "#", "_")
            // Truncate to 30 characters (Stata variable name limit minus "b_")
            if (strlen(nm) > 30) nm = substr(nm, 1, 30)
            // Prefix "b_" to ensure valid Stata name and avoid collision
            // with draw/sigma2/tau2
            all_names[3 + j] = "b_" + nm
        }
        else {
            all_names[3 + j] = "beta_" + strofreal(j)
        }
    }

    // --- Write to Stata dataset via preserve/restore ---
    // Strategy: preserve the current data, create the draws dataset, save it,
    // then restore the original data. This preserves the user's data in memory.
    //
    // IMPORTANT: This uses nested preserve/restore. If the calling code already
    // has a preserve active, this works correctly (Stata supports nested preserve).
    // However, do NOT call `horseshoe` with saving() from inside a Mata function
    // that has already called stata("preserve") — the restore here would pop the
    // wrong level. In normal usage (calling from a do-file), this is safe.
    stata("preserve")

    // Clear and set up the new dataset
    stata("clear")
    stata("quietly set obs " + strofreal(n_mcmc))

    // Create all variables in a single call (replaces 277 individual stata() calls)
    // st_addvar() requires a row vector, so transpose the column vector
    (void) st_addvar("double", all_names')

    // Store all data at once
    // st_store() takes (obs range, var indices, data matrix)
    st_store(., ., data)

    // Label the variables
    stata(`"label variable draw "MCMC draw number""')
    stata(`"label variable sigma2 "Posterior draw of sigma^2""')
    stata(`"label variable tau2 "Posterior draw of tau^2""')
    for (j = 1; j <= p; j++) {
        stata(`"label variable "' + all_names[3 + j] + ///
              `" "Posterior draw: "' + vnames[j] + `"""')
    }

    // Save the dataset (wrap filename in double quotes for paths with spaces)
    stata(`"quietly save ""' + filename + `"", replace"')

    // Restore original data
    stata("restore")
}


// -------------------------------------------------------------------------
// _hs_ess(chain)
//
// Compute effective sample size (ESS) for a single parameter chain using
// the initial positive sequence estimator (Geyer 1992). This is the same
// approach used by Stan's monitor() and R's coda::effectiveSize().
//
// Algorithm:
//   1. Compute autocorrelation at lags 0, 1, 2, ... via FFT-free method
//   2. Sum consecutive pairs of autocorrelations (lag 2k, lag 2k+1) until
//      the pair sum goes negative (initial positive sequence)
//   3. ESS = n / (1 + 2 * sum_of_autocorrelations)
//
// Arguments:
//   chain — n_mcmc × 1 column vector of draws
//
// Returns:
//   ESS as a real scalar (floored at 1 to avoid division-by-zero)
// -------------------------------------------------------------------------
real scalar _hs_ess(real colvector chain)
{
    real scalar    n, max_lag, k, rho_sum
    real scalar    rho_even, rho_odd, pair_sum
    real colvector centered
    real scalar    var0

    n = rows(chain)
    if (n < 4) return(n)

    // Center the chain
    centered = chain :- mean(chain)

    // Variance at lag 0 (denominator for autocorrelation)
    var0 = cross(centered, centered) / n
    if (var0 < 1e-300) return(n)  // constant chain — ESS = n

    // Compute autocorrelations and sum using initial positive sequence
    // We compute pairs: (rho_{2k}, rho_{2k+1}). Stop when pair sum < 0.
    max_lag = floor(n / 2) - 1
    rho_sum = 0

    for (k = 0; k <= max_lag; k++) {
        // Even lag: 2k
        rho_even = cross(centered[1::(n - 2*k)], centered[(2*k + 1)::n]) / (n * var0)

        // Odd lag: 2k + 1
        if (2*k + 1 < n) {
            rho_odd = cross(centered[1::(n - 2*k - 1)], centered[(2*k + 2)::n]) / (n * var0)
        }
        else {
            rho_odd = 0
        }

        // Initial positive sequence: stop if pair sum is negative
        pair_sum = rho_even + rho_odd
        if (pair_sum < 0) break

        rho_sum = rho_sum + pair_sum
    }

    // ESS = n / (1 + 2 * rho_sum)
    // The -1 adjusts for the lag-0 term already included in the pair at k=0
    // (rho_0 = 1 is always the first "even" term; we want 1 + 2*sum_{k>=1})
    // Actually: initial positive sequence sums pairs starting at (rho_0, rho_1),
    // so rho_sum already includes rho_0 = 1. The formula is:
    //   tau = -1 + 2 * rho_sum   (since rho_0 = 1 is counted, and we need
    //                              1 + 2*(rho_1 + rho_2 + ...) )
    // Equivalently: ESS = n / (-1 + 2 * rho_sum)
    // But more precisely: the initial positive sequence estimator gives
    //   tau_hat = -1 + 2 * sum_of_pair_sums
    // where pair_sum_k = rho_{2k} + rho_{2k+1}, and the first pair
    // includes rho_0 = 1. So tau_hat = -1 + 2 * rho_sum.
    //
    // ESS = n / tau_hat = n / (-1 + 2 * rho_sum)

    if (rho_sum < 1) rho_sum = 1  // guard: ensure ESS <= n

    return(max((1, n / (-1 + 2 * rho_sum))))
}


// -------------------------------------------------------------------------
// _hs_convergence_diag(res)
//
// Compute convergence diagnostics from posterior draws:
//   1. Effective sample size (ESS) for each beta_j, sigma2, tau2
//   2. Minimum ESS across all parameters
//   3. ESS for sigma2 and tau2 (key mixing indicators)
//
// Stores results as Stata scalars/matrices:
//   __hs_ess_beta   — 1 × p matrix: ESS for each beta_j
//   __hs_ess_min    — scalar: minimum ESS across all parameters
//   __hs_ess_sigma2 — scalar: ESS for sigma2
//   __hs_ess_tau2   — scalar: ESS for tau2
// -------------------------------------------------------------------------
void _hs_convergence_diag(struct _hs_result scalar res)
{
    real scalar    p, j, ess_j, ess_min
    real scalar    ess_sigma2, ess_tau2
    real rowvector ess_beta

    p = res.p

    // ESS for each beta_j
    ess_beta = J(1, p, .)
    ess_min = .

    for (j = 1; j <= p; j++) {
        ess_j = _hs_ess(res.beta_draws[j, .]')
        ess_beta[j] = ess_j
        if (ess_min == . | ess_j < ess_min) {
            ess_min = ess_j
        }
    }

    // ESS for sigma2 and tau2
    ess_sigma2 = _hs_ess(res.sigma2_draws)
    ess_tau2   = _hs_ess(res.tau2_draws)

    // Update minimum across all parameters
    if (ess_sigma2 < ess_min) ess_min = ess_sigma2
    if (ess_tau2   < ess_min) ess_min = ess_tau2

    // Store as Stata matrices/scalars
    st_matrix("__hs_ess_beta", ess_beta)
    st_numscalar("__hs_ess_min",    ess_min)
    st_numscalar("__hs_ess_sigma2", ess_sigma2)
    st_numscalar("__hs_ess_tau2",   ess_tau2)
}


// Note: _hs_quantile_pair() is provided by _hs_utils.ado (loaded at top)


// -------------------------------------------------------------------------
// _hs_stata_main(...)
//
// Main entry point called from the Stata ado-file.
// Orchestrates: data extraction → Gibbs sampler → post results → save draws
// -------------------------------------------------------------------------
void _hs_stata_main(
    string scalar depvar,
    string scalar indepvars,
    string scalar touse,
    string scalar pen_matname,
    real scalar   n_mcmc,
    real scalar   burnin,
    real scalar   thin,
    real scalar   lambda_scale,
    real scalar   tau_scale,
    real scalar   lower_pct,
    real scalar   upper_pct,
    real scalar   verbose,
    string scalar saving,
    string scalar varnames)
{
    // --- Extract data from Stata ---
    real colvector y
    real matrix    X
    real rowvector pen_flag
    real scalar    p

    // Get the estimation sample
    y = st_data(., depvar, touse)
    X = st_data(., tokens(indepvars), touse)
    p = cols(X)

    // --- Build penalization flag vector ---
    if (pen_matname != "") {
        // User supplied a penalized() matrix — read it
        pen_flag = st_matrix(pen_matname)
    }
    else {
        // Default: penalize all coefficients (standard horseshoe)
        pen_flag = J(1, p, 1)
    }

    // --- Run the Gibbs sampler ---
    struct _hs_result scalar res
    res = _hs_gibbs_run(y, X, pen_flag, n_mcmc, burnin, thin,
                        lambda_scale, tau_scale, verbose)

    // --- Post results to Stata e() ---
    _hs_post_results(res, lower_pct, upper_pct, lambda_scale, tau_scale,
                     burnin, thin, depvar, varnames)

    // --- Compute convergence diagnostics ---
    _hs_convergence_diag(res)

    // --- Save posterior draws if requested ---
    if (saving != "") {
        if (verbose) {
            printf("{txt}  [Horseshoe] Saving %g posterior draws to %s\n",
                   n_mcmc, saving)
            displayflush()
        }
        _hs_save_draws(res, saving, varnames)
    }

    if (verbose) {
        printf("{txt}  [Horseshoe Gibbs] Done.\n")
        displayflush()
    }
}


end
