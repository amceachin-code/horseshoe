# ===========================================================================
# R/hsreg_gibbs.R — Core Gibbs sampler for horseshoe regression
#
# Implements the auxiliary-variable Gibbs sampler from:
#   Makalic & Schmidt (2015), "A simple sampler for the horseshoe estimator"
#   arXiv:1508.03884v4
#
# Key feature: the `penalized` argument is a logical vector of length p.
#   penalized[j] = TRUE  --> beta_j gets the horseshoe prior (shrinkage)
#   penalized[j] = FALSE --> beta_j gets a flat (improper) prior (no shrinkage)
#
# All conditional posteriors are conjugate (see @details in roxygen).
#
# This is the raw sampler. For a full estimation workflow with S3 methods,
# see hsreg().
# ===========================================================================


#' Horseshoe Gibbs Sampler with Selective Penalization
#'
#' Samples from the posterior of Bayesian linear regression with the horseshoe
#' prior applied only to a user-specified subset of coefficients. Unpenalized
#' coefficients receive a flat (improper) prior — equivalent to an OLS-like
#' posterior for those parameters, conditional on the penalized parameters.
#'
#' This is the raw sampler that returns posterior draws directly. For a
#' higher-level interface with S3 methods (print, summary, predict, etc.),
#' see \code{\link{hsreg}}.
#'
#' @param y Numeric vector of length n. The outcome variable.
#' @param X Numeric matrix of dimension n x p. The design matrix. Include an
#'   intercept column explicitly if you want an unpenalized intercept.
#' @param penalized Logical vector of length p, or \code{NULL}. \code{TRUE}
#'   means the corresponding coefficient gets the horseshoe prior (shrinkage);
#'   \code{FALSE} means a flat prior (no shrinkage). Default (\code{NULL}):
#'   all coefficients are penalized (standard horseshoe).
#' @param lambda_scale Positive scalar. Scale of the half-Cauchy prior on the
#'   local shrinkage parameters: \eqn{\lambda_j \sim C^+(0, \text{lambda\_scale})}.
#'   Larger values allow individual coefficients to escape shrinkage more
#'   easily. Default: 1.
#' @param tau_scale Positive scalar. Scale of the half-Cauchy prior on the
#'   global shrinkage parameter: \eqn{\tau \sim C^+(0, \text{tau\_scale})}.
#'   Smaller values enforce more aggressive overall shrinkage. Default: 1.
#' @param n_mcmc Positive integer. Number of post-burnin MCMC draws to store.
#'   Default: 1000.
#' @param burnin Non-negative integer. Number of initial iterations to discard.
#'   Default: 500.
#' @param thin Positive integer. Keep every \code{thin}-th draw after burnin.
#'   Default: 1 (no thinning).
#' @param verbose Logical. Print progress every 500 iterations. Default: TRUE.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{beta_draws}{p x n_mcmc matrix of posterior draws for beta.}
#'     \item{sigma2_draws}{Numeric vector of length n_mcmc: posterior draws
#'       for \eqn{\sigma^2}.}
#'     \item{tau2_draws}{Numeric vector of length n_mcmc: posterior draws
#'       for \eqn{\tau^2}.}
#'     \item{lambda2_draws}{p_pen x n_mcmc matrix of posterior draws for
#'       \eqn{\lambda_j^2} (only for penalized coefficients, in the order
#'       they appear in X).}
#'     \item{pen_idx}{Integer vector: which columns of X are penalized.}
#'     \item{free_idx}{Integer vector: which columns of X are unpenalized.}
#'     \item{n_chol_fallback}{Integer: number of iterations where a ridge
#'       adjustment was needed for the Cholesky factorization.}
#'   }
#'
#' @details
#' The sampler follows Makalic & Schmidt (2015) Algorithm 1, with the
#' modification that the horseshoe hierarchy (local shrinkage \eqn{\lambda_j},
#' global shrinkage \eqn{\tau}, and their auxiliary variables \eqn{\nu_j},
#' \eqn{\xi}) only governs the penalized subset of coefficients. Unpenalized
#' coefficients enter the posterior of \eqn{\beta} through the likelihood
#' (\eqn{X^T X}) term only, with no prior penalty.
#'
#' The prior precision matrix \eqn{\Lambda^{*-1}} has:
#' \itemize{
#'   \item Entry \eqn{1/(\tau^2 \lambda_j^2)} at position (j,j) for penalized j
#'   \item Entry 0 at position (j,j) for unpenalized j (flat prior)
#' }
#'
#' The \eqn{\sigma^2} and \eqn{\tau^2} conditionals sum only over penalized
#' coefficients, since unpenalized betas contribute no prior information to
#' these hyperparameters.
#'
#' **Cholesky fallback**: If the Cholesky factorization of the precision
#' matrix \eqn{A = X^T X + \Lambda^{*-1}} fails (numerically singular), a
#' small ridge (\code{1e-8 * I_p}) is added and the factorization retried.
#' The count of such events is returned in \code{n_chol_fallback}. Persistent
#' fallbacks suggest near-collinear design columns.
#'
#' **Divergence detection**: After each sigma2 and tau2 draw, the sampler
#' checks for NaN or extreme values (> 1e200). If detected, the sampler
#' stops with an informative error message identifying the iteration.
#'
#' **Propriety requirement**: The posterior is proper when \eqn{n > p_{free}}
#' (number of unpenalized parameters) and the columns of X corresponding to
#' free parameters have full column rank.
#'
#' @references
#' Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe
#' estimator for sparse signals. \emph{Biometrika}, 97(2), 465-480.
#'
#' Makalic, E. & Schmidt, D. F. (2015). A simple sampler for the horseshoe
#' estimator. \emph{arXiv:1508.03884v4}.
#'
#' @examples
#' # Simple example: 5 true signals in 50 predictors
#' set.seed(123)
#' n <- 100; p <- 50
#' X <- matrix(rnorm(n * p), n, p)
#' beta_true <- c(rep(2, 5), rep(0, 45))
#' y <- X %*% beta_true + rnorm(n)
#'
#' fit <- hsreg_gibbs(y, X, n_mcmc = 200, burnin = 100, verbose = FALSE)
#' b_hat <- rowMeans(fit$beta_draws)
#' plot(beta_true, b_hat, xlab = "True", ylab = "Estimated")
#' abline(0, 1, col = "red")
#'
#' @export
hsreg_gibbs <- function(y, X, penalized = NULL,
                            lambda_scale = 1, tau_scale = 1,
                            n_mcmc = 1000, burnin = 500, thin = 1,
                            verbose = TRUE) {

  # -------------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------------
  if (!is.numeric(y)) stop("'y' must be a numeric vector")
  if (!is.matrix(X) || !is.numeric(X)) stop("'X' must be a numeric matrix")

  # Check for NA/NaN values — these propagate silently through crossprod
  if (anyNA(y)) stop("'y' contains NA values. Remove or impute before fitting.")
  if (anyNA(X)) stop("'X' contains NA values. Remove or impute before fitting.")

  n <- nrow(X)
  p <- ncol(X)

  if (length(y) != n) {
    stop("Length of y (", length(y), ") must equal nrow(X) (", n, ")")
  }

  # Default: penalize everything (standard horseshoe)
  if (is.null(penalized)) {
    penalized <- rep(TRUE, p)
  }
  if (!is.logical(penalized) || length(penalized) != p) {
    stop("'penalized' must be a logical vector of length ncol(X) (", p, ")")
  }

  # Indices for penalized vs free parameters
  pen_idx  <- which(penalized)
  free_idx <- which(!penalized)
  p_pen    <- length(pen_idx)
  p_free   <- length(free_idx)

  if (p_pen == 0L) {
    stop("No penalized parameters. Use OLS instead of the horseshoe sampler.")
  }

  # Check propriety: need n > p_free for the flat prior on free params
  if (p_free > 0L && n <= p_free) {
    warning("n (", n, ") <= p_free (", p_free,
            "). Posterior may be improper for unpenalized parameters. ",
            "Consider penalizing more columns or reducing p_free.")
  }

  # Validate scale parameters
  if (!is.numeric(lambda_scale) || length(lambda_scale) != 1L || lambda_scale <= 0) {
    stop("'lambda_scale' must be a positive scalar, got ", lambda_scale)
  }
  if (!is.numeric(tau_scale) || length(tau_scale) != 1L || tau_scale <= 0) {
    stop("'tau_scale' must be a positive scalar, got ", tau_scale)
  }

  # Validate MCMC parameters
  if (!is.numeric(n_mcmc) || length(n_mcmc) != 1L || n_mcmc < 1L) {
    stop("'n_mcmc' must be a positive integer, got ", n_mcmc)
  }
  if (!is.numeric(burnin) || length(burnin) != 1L || burnin < 0L) {
    stop("'burnin' must be a non-negative integer, got ", burnin)
  }
  if (!is.numeric(thin) || length(thin) != 1L || thin < 1L) {
    stop("'thin' must be a positive integer, got ", thin)
  }

  n_mcmc <- as.integer(n_mcmc)
  burnin <- as.integer(burnin)
  thin   <- as.integer(thin)

  # -------------------------------------------------------------------------
  # Precompute quantities that don't change across iterations
  # -------------------------------------------------------------------------
  XtX <- crossprod(X)        # p x p, symmetric positive semi-definite
  Xty <- crossprod(X, y)     # p x 1

  # Precompute 1/A^2 for the half-Cauchy scale-mixture auxiliary variables.
  # These appear in the conditional posteriors for nu_j and xi:
  #   nu_j | . ~ IG(1, 1/lambda_j^2 + 1/lambda_scale^2)
  #   xi   | . ~ IG(1, 1/tau^2     + 1/tau_scale^2)
  inv_lambda_scale_sq <- 1 / lambda_scale^2
  inv_tau_scale_sq    <- 1 / tau_scale^2

  # Total iterations: burnin + (n_mcmc draws * thinning interval)
  total_iter <- burnin + n_mcmc * thin

  # -------------------------------------------------------------------------
  # Initialize parameters
  # -------------------------------------------------------------------------
  # Beta: start at the ridge estimate (lambda = 1) for numerical stability.
  # This avoids the Cholesky failing at iteration 1 if X^T X is singular.
  beta <- as.numeric(solve(XtX + diag(1, p), Xty))

  # Residual variance: initialize from residuals
  sigma2 <- as.numeric(crossprod(y - X %*% beta)) / max(n - p, 1)

  # Local shrinkage: lambda_j^2 for penalized params only.
  # Initialize at 1 (no shrinkage). Length = p_pen.
  lambda2 <- rep(1, p_pen)

  # Global shrinkage: tau^2. Initialize at 1.
  tau2 <- 1

  # Auxiliary variables for the half-Cauchy decomposition:
  #   lambda_j ~ C+(0, lambda_scale) iff lambda_j^2 | nu_j ~ IG(1/2, 1/nu_j)
  #                                      and nu_j ~ IG(1/2, 1/lambda_scale^2)
  #   tau ~ C+(0, tau_scale)         iff tau^2 | xi ~ IG(1/2, 1/xi)
  #                                      and xi ~ IG(1/2, 1/tau_scale^2)
  # Initialize auxiliaries at 1.
  nu <- rep(1, p_pen)   # one per penalized param
  xi <- 1               # one for the global shrinkage

  # -------------------------------------------------------------------------
  # Storage for posterior draws
  # -------------------------------------------------------------------------
  beta_store    <- matrix(NA_real_, p, n_mcmc)
  sigma2_store  <- numeric(n_mcmc)
  tau2_store    <- numeric(n_mcmc)
  lambda2_store <- matrix(NA_real_, p_pen, n_mcmc)

  save_idx <- 0L  # counter for stored draws
  n_chol_fallback <- 0L  # counter for Cholesky ridge fallbacks

  # Identity matrix for potential ridge adjustment (allocated once)
  I_p <- diag(1, p)

  # -------------------------------------------------------------------------
  # Gibbs sampler main loop
  # -------------------------------------------------------------------------
  if (verbose) {
    message(sprintf("  [Horseshoe Gibbs] Starting sampler: %d total iterations (%d burnin + %d x %d thin)",
                    total_iter, burnin, n_mcmc, thin))
    message(sprintf("  [Horseshoe Gibbs] n = %d, p = %d (%d penalized, %d free)",
                    n, p, p_pen, p_free))
  }

  for (iter in seq_len(total_iter)) {

    # --- Step 1: Sample beta | rest  [eq. 9 of Makalic & Schmidt] ---
    #
    # beta | . ~ N_p(A^{-1} X^T y, sigma^2 * A^{-1})
    # where A = X^T X + Lambda_star_inv
    #
    # Lambda_star_inv is a p x p diagonal matrix:
    #   diagonal[j] = 1/(tau^2 * lambda_j^2)  for penalized j
    #   diagonal[j] = 0                        for unpenalized j (flat prior)
    diag_penalty <- numeric(p)
    diag_penalty[pen_idx] <- 1 / (tau2 * lambda2)
    A <- XtX + diag(diag_penalty)

    # Sample from N(mu, sigma^2 * A^{-1}) using Cholesky:
    #   A = R^T R  (R is upper triangular, from chol())
    #   mu = A^{-1} X^T y = solve(R, solve(t(R), Xty))
    #   beta = mu + sqrt(sigma2) * solve(R, z),  z ~ N(0, I_p)
    R <- tryCatch(chol(A), error = function(e) NULL)

    # Cholesky fallback: add small ridge and retry
    if (is.null(R)) {
      A <- A + 1e-8 * I_p
      R <- tryCatch(chol(A), error = function(e) NULL)
      n_chol_fallback <- n_chol_fallback + 1L
      if (is.null(R)) {
        stop("Cholesky factorization failed even with ridge adjustment at iteration ", iter,
             ". The design matrix may be severely rank-deficient.")
      }
    }

    mu_beta <- backsolve(R, forwardsolve(t(R), Xty))
    z <- stats::rnorm(p)
    beta <- as.numeric(mu_beta + sqrt(sigma2) * backsolve(R, z))

    # --- Step 2: Sample sigma^2 | rest  [eq. 10, modified] ---
    #
    # Only penalized betas contribute to the prior term.
    # sigma^2 | . ~ IG( (n + p_pen)/2,
    #                    [ ||y - X*beta||^2 + sum_{j in pen} beta_j^2/(tau^2*lambda_j^2) ] / 2 )
    resid <- y - X %*% beta
    rss <- as.numeric(crossprod(resid))
    prior_penalty <- sum(beta[pen_idx]^2 / (tau2 * lambda2))
    shape_sigma <- (n + p_pen) / 2
    rate_sigma  <- (rss + prior_penalty) / 2
    sigma2 <- 1 / stats::rgamma(1, shape = shape_sigma, rate = rate_sigma)

    # Divergence detection: check sigma2 for NaN/Inf
    if (is.nan(sigma2) || sigma2 > 1e200) {
      stop("Gibbs sampler diverged at iteration ", iter, ": sigma2 = ", sigma2,
           ". This usually indicates near-collinear design columns or ",
           "a highly ill-conditioned X'X matrix.")
    }

    # --- Step 3: Sample lambda_j^2 | rest  [eq. 11, penalized only] ---
    #
    # lambda_j^2 | . ~ IG(1, 1/nu_j + beta_j^2 / (2 * tau^2 * sigma^2))
    rate_lambda <- 1 / nu + beta[pen_idx]^2 / (2 * tau2 * sigma2)
    lambda2 <- 1 / stats::rgamma(p_pen, shape = 1, rate = rate_lambda)

    # --- Step 4: Sample tau^2 | rest  [eq. 12, summing over penalized only] ---
    #
    # tau^2 | . ~ IG( (p_pen + 1)/2,
    #                  1/xi + sum_{j in pen} beta_j^2/(lambda_j^2) / (2*sigma^2) )
    shape_tau <- (p_pen + 1) / 2
    rate_tau  <- 1 / xi + sum(beta[pen_idx]^2 / lambda2) / (2 * sigma2)
    tau2 <- 1 / stats::rgamma(1, shape = shape_tau, rate = rate_tau)

    # Divergence detection: check tau2 for NaN/Inf
    if (is.nan(tau2) || tau2 > 1e200) {
      stop("Gibbs sampler diverged at iteration ", iter, ": tau2 = ", tau2,
           ". This usually indicates a degenerate posterior for the ",
           "global shrinkage parameter.")
    }

    # --- Step 5: Sample auxiliary variables nu_j | rest  [penalized only] ---
    #
    # nu_j | . ~ IG(1, 1/lambda_scale^2 + 1/lambda_j^2)
    nu <- 1 / stats::rgamma(p_pen, shape = 1, rate = inv_lambda_scale_sq + 1 / lambda2)

    # --- Step 6: Sample auxiliary variable xi | rest ---
    #
    # xi | . ~ IG(1, 1/tau_scale^2 + 1/tau^2)
    xi <- 1 / stats::rgamma(1, shape = 1, rate = inv_tau_scale_sq + 1 / tau2)

    # --- Store draws (after burnin, respecting thinning) ---
    if (iter > burnin && (iter - burnin) %% thin == 0L) {
      save_idx <- save_idx + 1L
      beta_store[, save_idx]    <- beta
      sigma2_store[save_idx]    <- sigma2
      tau2_store[save_idx]      <- tau2
      lambda2_store[, save_idx] <- lambda2
    }

    # --- Progress reporting ---
    if (verbose && iter %% 500L == 0L) {
      message(sprintf("  [Horseshoe Gibbs] iter %d / %d | tau = %.4f | sigma = %.4f",
                      iter, total_iter, sqrt(tau2), sqrt(sigma2)))
    }
  }

  # --- Post-loop warnings ---
  if (n_chol_fallback > 0L && verbose) {
    warning("Cholesky factorization required ridge adjustment in ", n_chol_fallback,
            " of ", total_iter, " iterations. This may indicate near-collinear ",
            "columns in the design matrix.")
  }

  # -------------------------------------------------------------------------
  # Return posterior draws and index maps
  # -------------------------------------------------------------------------
  list(
    beta_draws      = beta_store,       # p x n_mcmc
    sigma2_draws    = sigma2_store,     # n_mcmc
    tau2_draws      = tau2_store,       # n_mcmc
    lambda2_draws   = lambda2_store,    # p_pen x n_mcmc
    pen_idx         = pen_idx,          # which columns were penalized
    free_idx        = free_idx,         # which columns were unpenalized
    n_chol_fallback = n_chol_fallback   # count of Cholesky ridge fallbacks
  )
}
