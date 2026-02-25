# ===========================================================================
# R/utils.R — Shared utilities for the horseshoe package
#
# Exports:
#   hs_ess()           — Effective sample size (initial positive sequence)
#   hs_quantile_pair() — Paired quantile in a single sort pass
#
# Internal:
#   .hs_convergence_diagnostics() — ESS for all parameters
# ===========================================================================


#' Effective Sample Size for an MCMC Chain
#'
#' Computes the effective sample size (ESS) of a single MCMC chain using the
#' initial positive sequence estimator (Geyer, 1992). This is the same approach
#' used by Stan's \code{monitor()} and R's \code{coda::effectiveSize()}.
#'
#' The ESS measures how many independent draws from the posterior an MCMC chain
#' of length \code{n} is equivalent to, accounting for autocorrelation between
#' successive draws. A chain with no autocorrelation has ESS = n; a highly
#' autocorrelated chain may have ESS << n.
#'
#' @param chain Numeric vector. A single MCMC chain (posterior draws for one
#'   parameter). Must have length >= 1.
#'
#' @return A scalar numeric value: the estimated effective sample size, floored
#'   at 1 (to avoid division-by-zero in downstream calculations). The ESS is
#'   always in the range \code{[1, length(chain)]}.
#'
#' @details
#' The algorithm works as follows:
#' \enumerate{
#'   \item Center the chain by subtracting its mean.
#'   \item Compute the lag-0 variance (denominator for autocorrelation).
#'   \item Sum consecutive pairs of autocorrelations:
#'     \eqn{(\rho_{2k}, \rho_{2k+1})} for \eqn{k = 0, 1, 2, \ldots}
#'   \item Stop when the pair sum goes negative (initial positive sequence
#'     criterion from Geyer, 1992).
#'   \item Return \eqn{n / (-1 + 2 \cdot \text{rho\_sum})}.
#' }
#'
#' For a constant chain (zero variance), ESS = n is returned, since every
#' draw is identical and there is no information loss from autocorrelation.
#'
#' @references
#' Geyer, C. J. (1992). Practical Markov chain Monte Carlo. \emph{Statistical
#' Science}, 7(4), 473-483.
#'
#' @examples
#' # IID draws — ESS should be close to n
#' set.seed(42)
#' iid_chain <- rnorm(1000)
#' hs_ess(iid_chain)
#'
#' # Highly autocorrelated chain — ESS << n
#' ar1_chain <- numeric(1000)
#' ar1_chain[1] <- rnorm(1)
#' for (i in 2:1000) ar1_chain[i] <- 0.99 * ar1_chain[i-1] + rnorm(1, sd = 0.1)
#' hs_ess(ar1_chain)
#'
#' @export
hs_ess <- function(chain) {

  # --- Input validation ---
  if (!is.numeric(chain)) {
    stop("'chain' must be a numeric vector")
  }

  n <- length(chain)

  # Short chains: return n directly (not enough data for autocorrelation)
  if (n < 4L) return(n)

  # --- Center the chain ---
  centered <- chain - mean(chain)

  # --- Variance at lag 0 (denominator for autocorrelation) ---
  var0 <- sum(centered^2) / n
  if (var0 < 1e-300) return(n)  # constant chain — ESS = n


  # --- Initial positive sequence estimator ---
  # Sum pairs of autocorrelations (rho_{2k}, rho_{2k+1}) until pair sum < 0.
  # This is Geyer's (1992) monotone positive sequence estimator truncated
  # at the first negative pair sum.
  max_lag <- floor(n / 2) - 1L
  rho_sum <- 0

  for (k in 0:max_lag) {
    # Even lag: 2k
    lag_even <- 2L * k
    idx_even_end <- n - lag_even
    rho_even <- sum(centered[1:idx_even_end] * centered[(lag_even + 1L):n]) / (n * var0)

    # Odd lag: 2k + 1
    lag_odd <- 2L * k + 1L
    if (lag_odd < n) {
      idx_odd_end <- n - lag_odd
      rho_odd <- sum(centered[1:idx_odd_end] * centered[(lag_odd + 1L):n]) / (n * var0)
    } else {
      rho_odd <- 0
    }

    # Initial positive sequence: stop if pair sum is negative
    pair_sum <- rho_even + rho_odd
    if (pair_sum < 0) break

    rho_sum <- rho_sum + pair_sum
  }

  # ESS = n / (-1 + 2 * rho_sum)
  # The -1 accounts for rho_0 = 1 being included in the first pair sum.
  # Guard: ensure ESS is in [1, n]
  if (rho_sum < 1) rho_sum <- 1
  ess <- n / (-1 + 2 * rho_sum)

  return(max(1, ess))
}


#' Paired Quantile Computation in a Single Sort Pass
#'
#' Computes two quantiles from the same numeric vector using a single sort
#' operation. This is more efficient than calling \code{\link[stats]{quantile}}
#' twice when you need both lower and upper bounds (e.g., for credible
#' intervals), because it avoids redundant sorting.
#'
#' @param v Numeric vector. The data from which to compute quantiles.
#' @param lo_prob Numeric scalar in \code{[0, 1]}. The lower quantile
#'   probability (e.g., 0.025 for a 95\% interval).
#' @param hi_prob Numeric scalar in \code{[0, 1]}. The upper quantile
#'   probability (e.g., 0.975 for a 95\% interval). Must be >= \code{lo_prob}.
#'
#' @return A named numeric vector of length 2: \code{c(lower = ..., upper = ...)}.
#'
#' @details
#' Uses R's default type-7 quantile interpolation (the same as
#' \code{stats::quantile(..., type = 7)}):
#' \deqn{h = (n - 1) \cdot p + 1}
#' where \eqn{p} is the probability and \eqn{n} is the length of the vector.
#' The quantile is linearly interpolated between the two bracketing sorted
#' values.
#'
#' @examples
#' set.seed(42)
#' x <- rnorm(1000)
#' hs_quantile_pair(x, 0.025, 0.975)
#'
#' # Verify it matches stats::quantile
#' all.equal(
#'   hs_quantile_pair(x, 0.025, 0.975),
#'   c(lower = unname(quantile(x, 0.025)),
#'     upper = unname(quantile(x, 0.975)))
#' )
#'
#' @export
hs_quantile_pair <- function(v, lo_prob, hi_prob) {

  # --- Input validation ---
  if (!is.numeric(v)) stop("'v' must be a numeric vector")
  n <- length(v)

  if (n == 0L) stop("'v' must have length >= 1")

  # Trivial case: single element

  if (n == 1L) {
    return(c(lower = v[1L], upper = v[1L]))
  }

  # --- Single sort ---
  sorted <- sort(v)

  # --- Lower quantile (R type-7 interpolation) ---
  h_lo <- (n - 1) * lo_prob + 1
  lo_lo <- max(1L, floor(h_lo))
  hi_lo <- min(n, ceiling(h_lo))
  if (lo_lo == hi_lo) {
    q_lower <- sorted[lo_lo]
  } else {
    q_lower <- sorted[lo_lo] + (h_lo - lo_lo) * (sorted[hi_lo] - sorted[lo_lo])
  }

  # --- Upper quantile (R type-7 interpolation) ---
  h_hi <- (n - 1) * hi_prob + 1
  lo_hi <- max(1L, floor(h_hi))
  hi_hi <- min(n, ceiling(h_hi))
  if (lo_hi == hi_hi) {
    q_upper <- sorted[lo_hi]
  } else {
    q_upper <- sorted[lo_hi] + (h_hi - lo_hi) * (sorted[hi_hi] - sorted[lo_hi])
  }

  return(c(lower = q_lower, upper = q_upper))
}


#' Convergence Diagnostics for Horseshoe Posterior Draws
#'
#' Computes effective sample size (ESS) for every parameter in the horseshoe
#' posterior: each beta_j, sigma2, and tau2. Returns a list of diagnostics
#' suitable for inclusion in a horseshoe fit object.
#'
#' @param beta_draws Numeric matrix, \code{p x n_mcmc}. Posterior draws for
#'   all regression coefficients.
#' @param sigma2_draws Numeric vector of length \code{n_mcmc}. Posterior draws
#'   for the residual variance.
#' @param tau2_draws Numeric vector of length \code{n_mcmc}. Posterior draws
#'   for the global shrinkage parameter.
#'
#' @return A list with components:
#'   \describe{
#'     \item{ess_beta}{Named numeric vector of length p: ESS for each beta_j.}
#'     \item{ess_sigma2}{Scalar: ESS for sigma2.}
#'     \item{ess_tau2}{Scalar: ESS for tau2.}
#'     \item{ess_min}{Scalar: minimum ESS across all parameters (beta, sigma2, tau2).}
#'   }
#'
#' @keywords internal
.hs_convergence_diagnostics <- function(beta_draws, sigma2_draws, tau2_draws) {

  p <- nrow(beta_draws)

  # ESS for each beta_j
  ess_beta <- apply(beta_draws, 1, hs_ess)

  # ESS for sigma2 and tau2
  ess_sigma2 <- hs_ess(sigma2_draws)
  ess_tau2 <- hs_ess(tau2_draws)

  # Minimum across all parameters
  ess_min <- min(c(ess_beta, ess_sigma2, ess_tau2))

  list(
    ess_beta   = ess_beta,
    ess_sigma2 = ess_sigma2,
    ess_tau2   = ess_tau2,
    ess_min    = ess_min
  )
}


#' Gelman-Rubin R-hat for a Single Parameter
#'
#' Computes the potential scale reduction factor (R-hat) from multiple
#' independent MCMC chains for a single scalar parameter. Values close to
#' 1.0 indicate convergence; values > 1.1 suggest the chains have not mixed.
#'
#' Uses the standard split-chain R-hat formula:
#' \deqn{W = \text{mean within-chain variance}}
#' \deqn{B = n \cdot \text{var(chain means)}}
#' \deqn{\hat{V} = (1 - 1/n) W + (1/n) B}
#' \deqn{\hat{R} = \sqrt{\hat{V} / W}}
#'
#' @param chain_list A list of numeric vectors, one per chain. Each vector
#'   contains the posterior draws for the same parameter from an independent
#'   chain. All chains must have the same length.
#'
#' @return A scalar numeric: the estimated R-hat. Returns \code{NaN} if
#'   within-chain variance is zero (constant chains).
#'
#' @keywords internal
.hs_rhat <- function(chain_list) {
  m <- length(chain_list)           # number of chains
  n <- length(chain_list[[1]])      # draws per chain

  # Within-chain variance: mean of per-chain variances
  chain_vars  <- vapply(chain_list, stats::var, numeric(1))
  W <- mean(chain_vars)

  # Between-chain variance: n * variance of chain means
  chain_means <- vapply(chain_list, mean, numeric(1))
  B <- n * stats::var(chain_means)

  # Pooled posterior variance estimate
  V_hat <- (1 - 1 / n) * W + (1 / n) * B

  # R-hat
  sqrt(V_hat / W)
}


#' Gelman-Rubin R-hat for Each Row of Per-Chain Draw Matrices
#'
#' Vectorized R-hat computation across all parameters (rows of the draw
#' matrices). Takes a list of \code{p x n_mcmc} matrices (one per chain)
#' and returns R-hat for each of the \code{p} parameters.
#'
#' @param per_chain_matrices A list of numeric matrices, each of dimension
#'   \code{p x n_mcmc}. One matrix per chain.
#'
#' @return A numeric vector of length \code{p}: R-hat for each parameter.
#'
#' @keywords internal
.hs_rhat_vector <- function(per_chain_matrices) {
  p <- nrow(per_chain_matrices[[1]])
  rhat <- numeric(p)

  for (j in seq_len(p)) {
    # Extract row j from each chain's matrix → list of numeric vectors
    chain_list_j <- lapply(per_chain_matrices, function(mat) mat[j, ])
    rhat[j] <- .hs_rhat(chain_list_j)
  }

  rhat
}
