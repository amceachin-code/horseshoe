*! _hs_utils.ado — Shared Mata utilities for the horseshoe package
*! Version 1.1.0 — 2026-02-21
*!
*! Provides shared Mata functions used by horseshoe.ado and horseshoe_cates.ado:
*!   _hs_quantile_pair()  — Paired lo/hi quantile in a single sort pass
*!
*! This file is loaded via `findfile _hs_utils.ado` + `run` from the
*! consuming .ado files. Placing shared code here eliminates duplication
*! and ensures bug fixes are applied in a single place.
*!
*! Author: Andrew McEachin
*! License: GPL-3.0

* Safely drop any existing definition before redefining
capture mata: mata drop _hs_quantile_pair()

mata:
mata set matastrict on


// -------------------------------------------------------------------------
// _hs_quantile_pair(v, lo_prob, hi_prob)
//
// Compute two quantiles from the same column vector with a single sort.
// Returns a 1×2 rowvector: (lo_quantile, hi_quantile).
//
// This is more efficient than calling a single-quantile function twice
// when you need both lower and upper bounds (avoids redundant sorting).
//
// Algorithm (R type 7 default):
//   h = (n - 1) * prob + 1   (1-based index into sorted vector)
//   lo = floor(h), hi = ceil(h)
//   quantile = sorted[lo] + (h - lo) * (sorted[hi] - sorted[lo])
// -------------------------------------------------------------------------
real rowvector _hs_quantile_pair(real colvector v,
                                 real scalar lo_prob,
                                 real scalar hi_prob)
{
    real colvector sorted
    real scalar    n, h_lo, h_hi, lo_lo, hi_lo, lo_hi, hi_hi
    real rowvector result

    n = rows(v)
    result = J(1, 2, .)

    if (n == 1) {
        result[1] = v[1]
        result[2] = v[1]
        return(result)
    }

    sorted = sort(v, 1)

    // Lower quantile
    h_lo = (n - 1) * lo_prob + 1
    lo_lo = floor(h_lo)
    hi_lo = ceil(h_lo)
    if (lo_lo < 1)  lo_lo = 1
    if (hi_lo > n)  hi_lo = n
    if (lo_lo > n)  lo_lo = n
    if (lo_lo == hi_lo) {
        result[1] = sorted[lo_lo]
    }
    else {
        result[1] = sorted[lo_lo] + (h_lo - lo_lo) * (sorted[hi_lo] - sorted[lo_lo])
    }

    // Upper quantile
    h_hi = (n - 1) * hi_prob + 1
    lo_hi = floor(h_hi)
    hi_hi = ceil(h_hi)
    if (lo_hi < 1)  lo_hi = 1
    if (hi_hi > n)  hi_hi = n
    if (lo_hi > n)  lo_hi = n
    if (lo_hi == hi_hi) {
        result[2] = sorted[lo_hi]
    }
    else {
        result[2] = sorted[lo_hi] + (h_hi - lo_hi) * (sorted[hi_hi] - sorted[lo_hi])
    }

    return(result)
}


end
