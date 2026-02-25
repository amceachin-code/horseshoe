# ===========================================================================
# R/hsreg.R — Main estimation function and S3 methods
#
# Provides:
#   hsreg()          — Fit Bayesian linear regression with horseshoe prior
#   print.hsreg()    — Display results table
#   summary.hsreg()  — Summary with coefficient table
#   coef.hsreg()     — Extract coefficient vector
#   vcov.hsreg()     — Extract posterior VCV matrix
#   predict.hsreg()  — Predictions (xb, residuals, stdp)
#   confint.hsreg()  — Credible intervals
# ===========================================================================


#' Fit Bayesian Linear Regression with the Horseshoe Prior
#'
#' Fits a Bayesian linear regression model using the horseshoe prior
#' (Carvalho, Polson & Scott, 2010) via the auxiliary-variable Gibbs sampler
#' of Makalic & Schmidt (2015). Supports selective per-parameter penalization:
#' user-specified coefficients can receive the horseshoe prior (shrinkage)
#' while others receive a flat (improper) prior (no shrinkage).
#'
#' Returns an S3 object of class \code{"hsreg"} with full method support:
#' \code{print}, \code{summary}, \code{coef}, \code{vcov}, \code{predict},
#' and \code{confint}.
#'
#' @param y Numeric vector of length n. The outcome variable.
#' @param X Numeric matrix of dimension n x p. The design matrix.
#' @param penalized Logical vector of length p, or \code{NULL}. \code{TRUE}
#'   = horseshoe prior (shrinkage); \code{FALSE} = flat prior. Default
#'   (\code{NULL}): all penalized.
#' @param lambda_scale Positive scalar. Scale of the half-Cauchy prior on
#'   local shrinkage. Default: 1.
#' @param tau_scale Positive scalar. Scale of the half-Cauchy prior on
#'   global shrinkage. Default: 1.
#' @param n_mcmc Positive integer. Number of post-burnin MCMC draws per chain.
#'   When \code{chains > 1}, total draws = \code{chains * n_mcmc}. Default: 1000.
#' @param burnin Non-negative integer. Burnin iterations (per chain). Default: 500.
#' @param thin Positive integer. Thinning interval (per chain). Default: 1.
#' @param level Numeric in (0, 1). Credible interval level. Default: 0.95.
#' @param scale_X Logical. If \code{TRUE}, standardize columns of X to unit
#'   variance before sampling, then rescale coefficients back to the original
#'   scale. Useful when columns of X have very different scales. Default: FALSE.
#' @param seed Integer or \code{NULL}. If not \code{NULL}, call
#'   \code{set.seed(seed)} before sampling for reproducibility. Default: NULL.
#' @param saving Character string or \code{NULL}. If not \code{NULL}, save
#'   the full posterior draws to a file. The format is detected from the
#'   extension: \code{.rds} (default) or \code{.csv}. Default: NULL.
#' @param varnames Character vector of length p, or \code{NULL}. Names for
#'   the coefficients. If \code{NULL}, uses \code{colnames(X)} or generates
#'   \code{X1, X2, ...}. Default: NULL.
#' @param verbose Logical. Print progress. Default: TRUE.
#' @param chains Positive integer. Number of independent MCMC chains. Each
#'   chain runs a full \code{n_mcmc} post-burnin draws. When \code{chains = 1}
#'   (default), behavior is bit-for-bit identical to single-chain mode. Per-chain
#'   seeds are derived from the master Mersenne Twister stream (not L'Ecuyer-CMRG),
#'   which is sufficient for typical use (2-10 chains). Default: 1.
#' @param cores Positive integer or \code{NULL}. Number of parallel workers for
#'   multi-chain runs. If \code{NULL} (default), auto-detects as
#'   \code{max(1, floor(detectCores(logical = FALSE) / 2))}, capped at
#'   \code{chains}. Uses PSOCK clusters for cross-platform compatibility
#'   (works on Windows, macOS, and Linux). When \code{cores = 1}, chains run
#'   sequentially with no parallel overhead. Default: NULL.
#'
#' @return An S3 object of class \code{"hsreg"}, which is a list with:
#'   \describe{
#'     \item{coefficients}{Named numeric vector of length p: posterior means.}
#'     \item{b_sd}{Named numeric vector: posterior standard deviations.}
#'     \item{b_lower}{Named numeric vector: lower credible interval bounds.}
#'     \item{b_upper}{Named numeric vector: upper credible interval bounds.}
#'     \item{vcov}{p x p posterior variance-covariance matrix.}
#'     \item{sigma2}{Scalar: posterior mean of \eqn{\sigma^2}.}
#'     \item{tau2}{Scalar: posterior mean of \eqn{\tau^2}.}
#'     \item{beta_draws}{p x n_mcmc matrix: full posterior draws for beta.}
#'     \item{sigma2_draws}{Numeric vector: posterior draws for \eqn{\sigma^2}.}
#'     \item{tau2_draws}{Numeric vector: posterior draws for \eqn{\tau^2}.}
#'     \item{lambda2_draws}{p_pen x n_mcmc matrix: posterior draws for \eqn{\lambda_j^2}.}
#'     \item{pen_idx, free_idx}{Integer vectors: penalized/free column indices.}
#'     \item{penalized}{Logical vector of length p: penalization flags.}
#'     \item{n, p, p_pen}{Integers: sample size, total predictors, penalized count.}
#'     \item{n_mcmc, burnin, thin, level}{MCMC settings.}
#'     \item{lambda_scale, tau_scale}{Prior scale parameters.}
#'     \item{X_sd}{Numeric vector or \code{NULL}: column SDs if \code{scale_X = TRUE}.}
#'     \item{varnames}{Character vector: coefficient names.}
#'     \item{seed, rng_state}{Seed and RNG state for reproducibility.}
#'     \item{n_chol_fallback}{Integer: Cholesky ridge fallback count.}
#'     \item{chains}{Integer: number of chains run.}
#'     \item{cores}{Integer: number of parallel workers used.}
#'     \item{total_draws}{Integer: total MCMC draws (\code{chains * n_mcmc}).}
#'     \item{rhat_beta}{Named numeric vector: Gelman-Rubin R-hat per coefficient.
#'       \code{NULL} if \code{chains = 1}.}
#'     \item{rhat_sigma2, rhat_tau2}{Scalar R-hat for variance parameters.
#'       \code{NULL} if \code{chains = 1}.}
#'     \item{rhat_max}{Scalar: maximum R-hat across all parameters.
#'       \code{NULL} if \code{chains = 1}.}
#'     \item{ess_beta}{Named numeric vector: ESS per coefficient.}
#'     \item{ess_sigma2, ess_tau2, ess_min}{ESS diagnostics.}
#'     \item{X, y}{Stored data (for \code{predict}).}
#'     \item{call}{The matched call.}
#'   }
#'
#' @references
#' Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe
#' estimator for sparse signals. \emph{Biometrika}, 97(2), 465-480.
#'
#' Makalic, E. & Schmidt, D. F. (2015). A simple sampler for the horseshoe
#' estimator. \emph{arXiv:1508.03884v4}.
#'
#' @examples
#' set.seed(42)
#' n <- 100; p <- 20
#' X <- matrix(rnorm(n * p), n, p)
#' beta_true <- c(rep(3, 5), rep(0, 15))
#' y <- X %*% beta_true + rnorm(n)
#'
#' fit <- hsreg(y, X, n_mcmc = 200, burnin = 100, verbose = FALSE)
#' print(fit)
#' predict(fit, type = "xb")[1:5]
#'
#' @export
hsreg <- function(y, X,
                      penalized = NULL,
                      lambda_scale = 1,
                      tau_scale = 1,
                      n_mcmc = 1000,
                      burnin = 500,
                      thin = 1,
                      level = 0.95,
                      scale_X = FALSE,
                      seed = NULL,
                      saving = NULL,
                      varnames = NULL,
                      verbose = TRUE,
                      chains = 1,
                      cores = NULL) {

  cl <- match.call()

  # -------------------------------------------------------------------------
  # Input validation (beyond what hsreg_gibbs checks)
  # -------------------------------------------------------------------------
  if (!is.numeric(y)) stop("'y' must be a numeric vector")
  if (!is.matrix(X) || !is.numeric(X)) stop("'X' must be a numeric matrix")

  n <- nrow(X)
  p <- ncol(X)

  if (length(y) != n) {
    stop("Length of y (", length(y), ") must equal nrow(X) (", n, ")")
  }

  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("'level' must be a scalar in (0, 1), got ", level)
  }

  # Default penalized: all TRUE
  if (is.null(penalized)) {
    penalized <- rep(TRUE, p)
  }
  if (!is.logical(penalized) || length(penalized) != p) {
    stop("'penalized' must be a logical vector of length ncol(X) (", p, ")")
  }

  # -------------------------------------------------------------------------
  # Validate chains / cores
  # -------------------------------------------------------------------------
  if (!is.numeric(chains) || length(chains) != 1L || chains < 1L) {
    stop("'chains' must be a positive integer, got ", chains)
  }
  chains <- as.integer(chains)

  if (is.null(cores)) {
    # Auto-detect: use half the physical cores, at least 1, capped at chains
    detected <- parallel::detectCores(logical = FALSE)
    if (is.na(detected)) detected <- 1L
    cores <- max(1L, floor(detected / 2))
    cores <- min(cores, chains)
  }
  if (!is.numeric(cores) || length(cores) != 1L || cores < 1L) {
    stop("'cores' must be a positive integer, got ", cores)
  }
  cores <- as.integer(cores)
  cores <- min(cores, chains)

  # -------------------------------------------------------------------------
  # Variable names
  # -------------------------------------------------------------------------
  if (is.null(varnames)) {
    if (!is.null(colnames(X))) {
      varnames <- colnames(X)
    } else {
      varnames <- paste0("X", seq_len(p))
    }
  }
  if (length(varnames) != p) {
    stop("'varnames' must have length ncol(X) (", p, "), got ", length(varnames))
  }

  # -------------------------------------------------------------------------
  # Seed handling — save and restore global RNG state so we don't
  # mutate the user's session. This follows the same pattern as
  # withr::with_seed().
  # -------------------------------------------------------------------------
  rng_state <- NULL
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = globalenv())
      } else {
        assign(".Random.seed", old_seed, envir = globalenv())
      }
    }, add = TRUE)
    set.seed(seed)
  }
  # Capture RNG state after set.seed (or current state if no seed)
  if (exists(".Random.seed", envir = globalenv())) {
    rng_state <- get(".Random.seed", envir = globalenv())
  }

  # -------------------------------------------------------------------------
  # Column scaling
  # -------------------------------------------------------------------------
  X_original <- X
  X_sd <- NULL

  if (scale_X) {
    X_sd <- apply(X, 2, stats::sd)
    # Replace zeros with 1 to avoid division by zero (constant columns)
    X_sd[X_sd == 0] <- 1
    X <- sweep(X, 2, X_sd, FUN = "/")
  }

  # -------------------------------------------------------------------------
  # Run the Gibbs sampler (single-chain or multi-chain)
  # -------------------------------------------------------------------------

  if (chains == 1L) {
    # ----- Single chain: bit-for-bit identical to original code path -----
    fit <- hsreg_gibbs(
      y = y, X = X, penalized = penalized,
      lambda_scale = lambda_scale, tau_scale = tau_scale,
      n_mcmc = n_mcmc, burnin = burnin, thin = thin,
      verbose = verbose
    )

    beta_draws    <- fit$beta_draws       # p x n_mcmc
    sigma2_draws  <- fit$sigma2_draws     # n_mcmc vector
    tau2_draws    <- fit$tau2_draws       # n_mcmc vector
    lambda2_draws <- fit$lambda2_draws    # p_pen x n_mcmc
    pen_idx       <- fit$pen_idx
    free_idx      <- fit$free_idx
    n_chol_fallback <- fit$n_chol_fallback

    # R-hat: not computable from a single chain
    rhat_beta   <- NULL
    rhat_sigma2 <- NULL
    rhat_tau2   <- NULL
    rhat_max    <- NULL

  } else {
    # ----- Multi-chain: generate per-chain seeds from master RNG -----
    chain_seeds <- sample.int(.Machine$integer.max, chains)

    # Capture the hsreg_gibbs function object for export to workers.
    # We export the function object directly rather than using
    # loadNamespace("hsreg") so this works during devtools::load_all().
    gibbs_fn <- hsreg_gibbs

    # Worker function: runs one chain with its own seed
    run_one_chain <- function(chain_id) {
      set.seed(chain_seeds[chain_id])
      gibbs_fn(
        y = y, X = X, penalized = penalized,
        lambda_scale = lambda_scale, tau_scale = tau_scale,
        n_mcmc = n_mcmc, burnin = burnin, thin = thin,
        verbose = FALSE
      )
    }

    if (cores > 1L) {
      # --- Parallel dispatch via PSOCK cluster (cross-platform) ---
      if (verbose) {
        message(sprintf("  [Horseshoe] Launching %d chains on %d cores...",
                        chains, cores))
      }
      cl_par <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl_par), add = TRUE)

      # Export all objects needed by the worker function
      parallel::clusterExport(cl_par,
        varlist = c("chain_seeds", "gibbs_fn", "y", "X", "penalized",
                     "lambda_scale", "tau_scale", "n_mcmc", "burnin", "thin"),
        envir = environment()
      )

      chain_fits <- parallel::parLapply(cl_par, seq_len(chains), run_one_chain)

    } else {
      # --- Sequential multi-chain (cores == 1) ---
      chain_fits <- vector("list", chains)
      for (cc in seq_len(chains)) {
        if (verbose) {
          message(sprintf("  [Horseshoe] Running chain %d / %d...", cc, chains))
        }
        chain_fits[[cc]] <- run_one_chain(cc)
      }
    }

    # ----- R-hat diagnostics (computed BEFORE combining to avoid double memory) -----
    rhat_beta <- .hs_rhat_vector(lapply(chain_fits, "[[", "beta_draws"))
    names(rhat_beta) <- varnames

    rhat_sigma2 <- .hs_rhat(lapply(chain_fits, "[[", "sigma2_draws"))
    rhat_tau2   <- .hs_rhat(lapply(chain_fits, "[[", "tau2_draws"))
    rhat_max    <- max(c(rhat_beta, rhat_sigma2, rhat_tau2), na.rm = TRUE)

    if (rhat_max > 1.1 && verbose) {
      warning("Max R-hat = ", round(rhat_max, 3),
              " (> 1.1). Chains may not have converged. ",
              "Consider increasing n_mcmc or burnin.")
    }

    # ----- Combine draws across chains -----
    beta_draws    <- do.call(cbind, lapply(chain_fits, "[[", "beta_draws"))
    sigma2_draws  <- unlist(lapply(chain_fits, "[[", "sigma2_draws"))
    tau2_draws    <- unlist(lapply(chain_fits, "[[", "tau2_draws"))
    lambda2_draws <- do.call(cbind, lapply(chain_fits, "[[", "lambda2_draws"))

    # Metadata from first chain (all chains share the same structure)
    pen_idx   <- chain_fits[[1]]$pen_idx
    free_idx  <- chain_fits[[1]]$free_idx
    n_chol_fallback <- sum(vapply(chain_fits, "[[", numeric(1), "n_chol_fallback"))
  }

  total_draws <- ncol(beta_draws)

  # -------------------------------------------------------------------------
  # Post-process: compute posterior summaries (works on combined draws)
  # -------------------------------------------------------------------------

  # If scaled: rescale beta draws back to original parameterization
  # beta_orig = beta_scaled / X_sd
  if (scale_X) {
    beta_draws <- beta_draws / X_sd  # vectorized: each row divided by X_sd[j]
  }

  # Posterior means
  b_mean <- rowMeans(beta_draws)

  # Posterior SDs
  b_sd <- apply(beta_draws, 1, stats::sd)

  # Credible intervals from quantiles
  alpha <- 1 - level
  lo_prob <- alpha / 2
  hi_prob <- 1 - alpha / 2

  b_lower <- numeric(p)
  b_upper <- numeric(p)
  for (j in seq_len(p)) {
    qi <- hs_quantile_pair(beta_draws[j, ], lo_prob, hi_prob)
    b_lower[j] <- qi["lower"]
    b_upper[j] <- qi["upper"]
  }

  # Name the vectors
  names(b_mean)  <- varnames
  names(b_sd)    <- varnames
  names(b_lower) <- varnames
  names(b_upper) <- varnames

  # Posterior variance-covariance matrix
  # V = Var(beta_draws^T) where beta_draws is p x total_draws
  # stats::var expects observations in rows, variables in columns
  V <- stats::var(t(beta_draws))
  rownames(V) <- varnames
  colnames(V) <- varnames

  # -------------------------------------------------------------------------
  # Convergence diagnostics (ESS on combined draws)
  # -------------------------------------------------------------------------
  diag <- .hs_convergence_diagnostics(beta_draws, sigma2_draws, tau2_draws)
  names(diag$ess_beta) <- varnames

  # Warn if ESS is suspiciously low
  if (diag$ess_min < 100 && verbose) {
    warning("Minimum ESS = ", round(diag$ess_min),
            " (< 100). Consider increasing n_mcmc or burnin.")
  }

  # -------------------------------------------------------------------------
  # Save draws if requested
  # -------------------------------------------------------------------------
  if (!is.null(saving)) {
    draws_df <- data.frame(
      draw = seq_len(total_draws),
      sigma2 = sigma2_draws,
      tau2 = tau2_draws,
      t(beta_draws)
    )
    colnames(draws_df) <- c("draw", "sigma2", "tau2", varnames)

    ext <- tolower(tools::file_ext(saving))
    if (ext == "csv") {
      utils::write.csv(draws_df, file = saving, row.names = FALSE)
    } else {
      # Default to .rds
      saveRDS(draws_df, file = saving)
    }
    if (verbose) {
      message("  [Horseshoe] Posterior draws saved to: ", saving)
    }
  }

  # -------------------------------------------------------------------------
  # Build and return S3 object
  # -------------------------------------------------------------------------
  result <- structure(list(
    # Point estimates and summaries
    coefficients  = b_mean,
    b_sd          = b_sd,
    b_lower       = b_lower,
    b_upper       = b_upper,
    vcov          = V,
    sigma2        = mean(sigma2_draws),
    tau2          = mean(tau2_draws),

    # Full posterior draws (for custom analysis)
    beta_draws    = beta_draws,        # p x total_draws (rescaled if scale_X)
    sigma2_draws  = sigma2_draws,
    tau2_draws    = tau2_draws,
    lambda2_draws = lambda2_draws,

    # Penalization info
    pen_idx       = pen_idx,
    free_idx      = free_idx,
    penalized     = penalized,

    # Metadata
    n = n, p = p, p_pen = length(pen_idx),
    n_mcmc = n_mcmc, burnin = burnin, thin = thin,
    level = level,
    lambda_scale = lambda_scale, tau_scale = tau_scale,
    X_sd = X_sd, varnames = varnames,
    seed = seed, rng_state = rng_state,
    n_chol_fallback = n_chol_fallback,

    # Multi-chain metadata
    chains      = chains,
    cores       = cores,
    total_draws = total_draws,

    # R-hat diagnostics (NULL for single chain)
    rhat_beta   = rhat_beta,
    rhat_sigma2 = rhat_sigma2,
    rhat_tau2   = rhat_tau2,
    rhat_max    = rhat_max,

    # Convergence diagnostics (ESS on combined draws)
    ess_beta   = diag$ess_beta,
    ess_sigma2 = diag$ess_sigma2,
    ess_tau2   = diag$ess_tau2,
    ess_min    = diag$ess_min,

    # Stored data (for predict)
    X = X_original, y = y,
    call = cl
  ), class = "hsreg")

  return(result)
}


# ===========================================================================
# S3 Methods
# ===========================================================================


#' Print Method for Horseshoe Objects
#'
#' Displays a formatted results table mirroring the Stata \code{hsreg}
#' output, including model summary statistics, convergence diagnostics, and
#' a coefficient table with posterior means, SDs, credible intervals, and
#' penalization status.
#'
#' @param x An object of class \code{"hsreg"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.hsreg <- function(x, ...) {

  # --- Header block ---
  cat("\n")
  cat(format("Horseshoe prior regression", width = 55),
      sprintf("Number of obs   = %8d\n", x$n))
  cat(format("", width = 55),
      sprintf("Predictors      = %8d\n", x$p))
  cat(format("", width = 55),
      sprintf("Penalized       = %8d\n", x$p_pen))
  if (x$chains > 1L) {
    cat(format("", width = 55),
        sprintf("Chains          = %8d\n", x$chains))
    cat(format("", width = 55),
        sprintf("Draws/chain     = %8d\n", x$n_mcmc))
    cat(format("", width = 55),
        sprintf("Total draws     = %8d\n", x$total_draws))
  } else {
    cat(format("", width = 55),
        sprintf("MCMC draws      = %8d\n", x$n_mcmc))
  }
  cat("\n")

  # --- Hyperparameter summaries ---
  cat(sprintf("Posterior mean of sigma^2 = %9.4f\n", x$sigma2))
  cat(sprintf("Posterior mean of tau^2   = %9.6f\n", x$tau2))
  cat("\n")

  # --- Convergence diagnostics ---
  cat(sprintf("Min effective sample size = %9.0f\n", x$ess_min))
  cat(sprintf("ESS for sigma^2           = %9.0f\n", x$ess_sigma2))
  cat(sprintf("ESS for tau^2             = %9.0f\n", x$ess_tau2))

  if (!is.null(x$rhat_max)) {
    cat(sprintf("Max R-hat                 = %9.3f\n", x$rhat_max))
    if (x$rhat_max > 1.1) {
      cat("Warning: max R-hat > 1.1. Chains may not have converged.\n")
    }
  }

  if (x$ess_min < 100) {
    cat("Warning: minimum ESS < 100. Consider increasing n_mcmc or burnin.\n")
  }
  cat("\n")

  # --- Coefficient table ---
  # Header
  sep <- paste0(strrep("-", 14), "|", strrep("-", 61))
  cat(sep, "\n")
  cat(sprintf("%13s | %11s %11s %11s %11s %10s\n",
              "Variable", "Post.Mean", "Post.SD", "Lower CI", "Upper CI", "Penalized"))
  cat(sep, "\n")

  # Body — one row per coefficient
  for (j in seq_len(x$p)) {
    # Truncate long variable names
    vname <- substr(x$varnames[j], 1, 12)
    pen_str <- if (x$penalized[j]) "Yes" else "No"
    cat(sprintf("%13s | %11.4f %11.4f %11.4f %11.4f %10s\n",
                vname,
                x$coefficients[j],
                x$b_sd[j],
                x$b_lower[j],
                x$b_upper[j],
                pen_str))
  }

  cat(sep, "\n")
  cat(sprintf("CI level: %g%%\n", x$level * 100))

  invisible(x)
}


#' Summary Method for Horseshoe Objects
#'
#' Returns a summary object containing a data.frame of coefficient estimates
#' with posterior means, standard deviations, credible intervals, and
#' penalization flags. The summary object has its own print method.
#'
#' @param object An object of class \code{"hsreg"}.
#' @param level Numeric in (0, 1). Credible interval level. If different from
#'   the level used at estimation, intervals are recomputed from the stored
#'   posterior draws. Default: uses the level from the fit.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class \code{"summary.hsreg"} containing:
#'   \describe{
#'     \item{coef_table}{A data.frame with columns: Post.Mean, Post.SD,
#'       Lower, Upper, Penalized.}
#'     \item{sigma2}{Posterior mean of sigma^2.}
#'     \item{tau2}{Posterior mean of tau^2.}
#'     \item{n, p, p_pen, n_mcmc}{Model dimensions and MCMC settings.}
#'     \item{ess_min, ess_sigma2, ess_tau2}{Convergence diagnostics.}
#'     \item{level}{The CI level used.}
#'   }
#'
#' @export
summary.hsreg <- function(object, level = object$level, ...) {

  # Recompute CIs if level differs from original

  if (abs(level - object$level) > 1e-10) {
    alpha <- 1 - level
    lo_prob <- alpha / 2
    hi_prob <- 1 - alpha / 2
    b_lower <- numeric(object$p)
    b_upper <- numeric(object$p)
    for (j in seq_len(object$p)) {
      qi <- hs_quantile_pair(object$beta_draws[j, ], lo_prob, hi_prob)
      b_lower[j] <- qi["lower"]
      b_upper[j] <- qi["upper"]
    }
  } else {
    b_lower <- object$b_lower
    b_upper <- object$b_upper
  }

  coef_table <- data.frame(
    Post.Mean = unname(object$coefficients),
    Post.SD   = unname(object$b_sd),
    Lower     = unname(b_lower),
    Upper     = unname(b_upper),
    Penalized = object$penalized,
    row.names = object$varnames,
    stringsAsFactors = FALSE
  )

  result <- structure(list(
    coef_table   = coef_table,
    sigma2       = object$sigma2,
    tau2         = object$tau2,
    n            = object$n,
    p            = object$p,
    p_pen        = object$p_pen,
    n_mcmc       = object$n_mcmc,
    chains       = object$chains,
    total_draws  = object$total_draws,
    rhat_max     = object$rhat_max,
    ess_min      = object$ess_min,
    ess_sigma2   = object$ess_sigma2,
    ess_tau2     = object$ess_tau2,
    level        = level
  ), class = "summary.hsreg")

  return(result)
}


#' Print Method for Summary of Horseshoe Objects
#'
#' @param x An object of class \code{"summary.hsreg"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.summary.hsreg <- function(x, ...) {
  cat("\nHorseshoe prior regression summary\n\n")
  if (!is.null(x$chains) && x$chains > 1L) {
    cat(sprintf("  n = %d, p = %d (%d penalized), %d chains x %d = %d total draws\n",
                x$n, x$p, x$p_pen, x$chains, x$n_mcmc, x$total_draws))
  } else {
    cat(sprintf("  n = %d, p = %d (%d penalized), n_mcmc = %d\n",
                x$n, x$p, x$p_pen, x$n_mcmc))
  }
  cat(sprintf("  sigma^2 = %.4f, tau^2 = %.6f\n", x$sigma2, x$tau2))
  cat(sprintf("  Min ESS = %.0f, ESS(sigma^2) = %.0f, ESS(tau^2) = %.0f\n",
              x$ess_min, x$ess_sigma2, x$ess_tau2))
  if (!is.null(x$rhat_max)) {
    cat(sprintf("  Max R-hat = %.3f\n", x$rhat_max))
  }
  cat("\n")
  cat(sprintf("Coefficients (%g%% credible intervals):\n", x$level * 100))
  print(x$coef_table)
  invisible(x)
}


#' Extract Coefficients from a Horseshoe Object
#'
#' @param object An object of class \code{"hsreg"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Named numeric vector of posterior mean coefficients.
#'
#' @export
coef.hsreg <- function(object, ...) {
  object$coefficients
}


#' Extract Posterior Variance-Covariance Matrix
#'
#' @param object An object of class \code{"hsreg"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A p x p numeric matrix: the posterior variance-covariance matrix
#'   of the regression coefficients (computed from the MCMC draws).
#'
#' @export
vcov.hsreg <- function(object, ...) {
  object$vcov
}


#' Predictions from a Horseshoe Object
#'
#' Compute fitted values, residuals, or standard errors of the linear
#' predictor from a fitted horseshoe model.
#'
#' @param object An object of class \code{"hsreg"}.
#' @param newdata Numeric matrix with p columns. New design matrix for
#'   prediction. If \code{NULL}, uses the stored training data. Default: NULL.
#' @param type Character string specifying the type of prediction:
#'   \describe{
#'     \item{\code{"xb"}}{Linear predictor: \code{newdata \%*\% coef(object)}.
#'       This is the default.}
#'     \item{\code{"residuals"}}{Residuals: \code{y - fitted}. Requires
#'       \code{newdata = NULL} (uses stored training data).}
#'     \item{\code{"stdp"}}{Standard error of the linear predictor:
#'       \code{sqrt(diag(newdata \%*\% V \%*\% t(newdata)))}.}
#'   }
#' @param ... Additional arguments (currently ignored).
#'
#' @return Numeric vector of predictions.
#'
#' @export
predict.hsreg <- function(object, newdata = NULL,
                               type = c("xb", "residuals", "stdp"), ...) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    X_pred <- object$X
  } else {
    if (!is.matrix(newdata) || !is.numeric(newdata)) {
      stop("'newdata' must be a numeric matrix")
    }
    if (ncol(newdata) != object$p) {
      stop("'newdata' must have ", object$p, " columns, got ", ncol(newdata))
    }
    X_pred <- newdata
  }

  b <- object$coefficients

  if (type == "xb") {
    return(as.numeric(X_pred %*% b))
  }

  if (type == "residuals") {
    if (!is.null(newdata)) {
      stop("type = 'residuals' requires newdata = NULL (uses stored training data)")
    }
    fitted <- as.numeric(X_pred %*% b)
    return(object$y - fitted)
  }

  if (type == "stdp") {
    # Standard error of the linear predictor
    # stdp_i = sqrt( x_i' V x_i ) for each row x_i
    V <- object$vcov
    # Efficient computation: (X %*% V) * X, then rowSums
    XVt <- X_pred %*% V
    stdp <- sqrt(pmax(0, rowSums(XVt * X_pred)))
    return(stdp)
  }
}


#' Credible Intervals for Horseshoe Coefficients
#'
#' Compute Bayesian credible intervals for the regression coefficients
#' from a fitted horseshoe model. If the requested level differs from the
#' level used at estimation, intervals are recomputed from the stored
#' posterior draws.
#'
#' @param object An object of class \code{"hsreg"}.
#' @param parm Character or integer vector specifying which parameters to
#'   return intervals for. If \code{NULL}, all parameters are returned.
#'   Default: NULL.
#' @param level Numeric in (0, 1). Credible interval level. Default: uses
#'   the level from the fit.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A matrix with columns for the lower and upper bounds, and rows
#'   named by the coefficient names.
#'
#' @export
confint.hsreg <- function(object, parm = NULL, level = object$level, ...) {

  # Recompute if level differs
  if (abs(level - object$level) > 1e-10) {
    alpha <- 1 - level
    lo_prob <- alpha / 2
    hi_prob <- 1 - alpha / 2
    b_lower <- numeric(object$p)
    b_upper <- numeric(object$p)
    for (j in seq_len(object$p)) {
      qi <- hs_quantile_pair(object$beta_draws[j, ], lo_prob, hi_prob)
      b_lower[j] <- qi["lower"]
      b_upper[j] <- qi["upper"]
    }
    names(b_lower) <- object$varnames
    names(b_upper) <- object$varnames
  } else {
    b_lower <- object$b_lower
    b_upper <- object$b_upper
  }

  # Build CI matrix
  pct_lo <- sprintf("%.1f %%", (1 - level) / 2 * 100)
  pct_hi <- sprintf("%.1f %%", (1 - (1 - level) / 2) * 100)

  ci <- cbind(b_lower, b_upper)
  colnames(ci) <- c(pct_lo, pct_hi)
  rownames(ci) <- object$varnames

  # Subset if parm is specified

  if (!is.null(parm)) {
    if (is.character(parm)) {
      idx <- match(parm, object$varnames)
      if (any(is.na(idx))) {
        stop("Unknown parameter names: ",
             paste(parm[is.na(idx)], collapse = ", "))
      }
    } else {
      idx <- parm
    }
    ci <- ci[idx, , drop = FALSE]
  }

  return(ci)
}
