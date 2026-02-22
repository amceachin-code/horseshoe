*! horseshoe_p.ado — Predict after horseshoe
*! Version 1.1.0 — 2026-02-21
*!
*! Postestimation predict command for horseshoe.ado.
*! Supports: xb (default), residuals, stdp (via e(V)).
*!
*! This file follows the pattern of rreg_p.ado (robust regression predict).
*! It is registered via  ereturn local predict "horseshoe_p"
*! so that  predict newvar [, options]  dispatches here automatically.
*!
*! Author: Andrew McEachin
*! License: GPL-3.0

program define horseshoe_p
    version 16.0

    * -------------------------------------------------------------------
    * SYNTAX PARSING
    *
    * Uses Stata's standard syntax command. Option names:
    *   XB        — linear prediction (default)
    *   Residuals — y - xb
    *   STDP      — standard error of prediction (requires e(V))
    *
    * The capital letters define the minimum abbreviation:
    *   xb, r[esiduals], stdp
    * -------------------------------------------------------------------
    syntax newvarname [if] [in] [, XB Residuals STDP]

    * Verify this is being called after horseshoe
    if "`e(cmd)'" != "horseshoe" {
        display as error "horseshoe_p requires prior estimation by {bf:horseshoe}"
        error 301
    }

    * Count how many options were specified (at most one allowed)
    local nopt : word count `xb' `residuals' `stdp'
    if `nopt' > 1 {
        display as error "only one option allowed"
        exit 198
    }

    * -------------------------------------------------------------------
    * STDP — standard error of the linear prediction
    *
    * Delegates to _predict which computes sqrt(x' * V * x) using e(V).
    * -------------------------------------------------------------------
    if "`stdp'" != "" {
        _predict `typlist' `varlist' `if' `in', stdp
        label variable `varlist' "S.E. of prediction"
        exit
    }

    * -------------------------------------------------------------------
    * RESIDUALS — observed minus fitted: y - xb
    * -------------------------------------------------------------------
    if "`residuals'" != "" {
        tempvar xb_tmp
        quietly _predict double `xb_tmp' `if' `in'
        generate `typlist' `varlist' = `e(depvar)' - `xb_tmp' `if' `in'
        label variable `varlist' "Residuals"
        exit
    }

    * -------------------------------------------------------------------
    * XB (default) — linear prediction X*b
    * -------------------------------------------------------------------
    _predict `typlist' `varlist' `if' `in'
    label variable `varlist' "Fitted values"
end
