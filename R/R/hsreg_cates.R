# ===========================================================================
# R/hsreg_cates.R — CATE extraction for multi-arm treatment designs
#
# Provides:
#   hsreg_cates()           — Extract CATEs from a fitted horseshoe object
#   hsreg_cates_from_file() — Extract CATEs from a saved draws file
#
# These functions are general-purpose tools for extracting individual-level
# conditional average treatment effects (iCATEs) from posterior beta draws
# when the design matrix follows the saturated interaction layout:
#   [X(n_x), D(n_d), X:D(n_d * n_x)]
#
# No ACIC-specific or project-specific code is included.
# ===========================================================================


#' Extract CATEs from a Fitted Horseshoe Object
#'
#' Extracts individual-level conditional average treatment effects (iCATEs)
#' from posterior draws of a horseshoe regression with a multi-arm treatment
#' design. The design matrix must follow the saturated interaction layout:
#' \code{[X(n_x), D(n_d), X:D(n_d * n_x)]}.
#'
#' For each treatment contrast \eqn{w = 1, \ldots, n_d}, the iCATE for
#' observation \eqn{i} is:
#' \deqn{\tau_w(x_i) = \gamma_w + x_i' \delta_w}
#' where \eqn{\gamma_w} is the treatment dummy coefficient and
#' \eqn{\delta_w} is the vector of treatment-covariate interaction
#' coefficients.
#'
#' @param fit An object of class \code{"hsreg"}, or a list with a
#'   \code{beta_draws} element (p x n_mcmc matrix).
#' @param X_test Numeric matrix of dimension n_test x n_x. Test covariate
#'   values at which to evaluate the CATEs. These should be the raw
#'   (unscaled) covariates.
#' @param n_x Integer. Number of covariate columns in the design matrix.
#' @param n_d Integer. Number of treatment dummy columns (contrasts).
#' @param X_sd Numeric vector of length p, or \code{NULL}. Column standard
#'   deviations used for scaling. If \code{NULL} and \code{fit} has an
#'   \code{X_sd} element, that is used. If no scaling was applied, leave
#'   as \code{NULL}. Default: NULL.
#' @param level Numeric in (0, 1). Credible interval level. Default: 0.95.
#'
#' @return A list with components:
#'   \describe{
#'     \item{cate_hat}{n_test x n_d matrix: posterior mean iCATEs.}
#'     \item{cate_lo}{n_test x n_d matrix: lower credible interval bounds.}
#'     \item{cate_hi}{n_test x n_d matrix: upper credible interval bounds.}
#'     \item{n_test}{Integer: number of test observations.}
#'     \item{n_d}{Integer: number of treatment contrasts.}
#'     \item{n_draws}{Integer: number of MCMC draws used.}
#'     \item{level}{Numeric: CI level used.}
#'   }
#'
#' @details
#' **Design matrix layout assumption**: The function assumes the design
#' matrix was constructed with columns in this order:
#' \enumerate{
#'   \item Columns 1 to n_x: Covariate main effects (X)
#'   \item Columns n_x+1 to n_x+n_d: Treatment dummies (D)
#'   \item Columns n_x+n_d+1 to end: Interactions (X:D), blocked by arm.
#'     That is, n_x columns for contrast 1, then n_x for contrast 2, etc.
#' }
#'
#' The total number of columns must equal \code{n_x + n_d + n_d * n_x}.
#'
#' **Rescaling**: If the horseshoe model was fit on a scaled design matrix
#' (\code{scale_X = TRUE}), the beta draws are in the scaled parameterization.
#' This function rescales them back to the original scale using \code{X_sd}
#' before computing CATEs. If \code{fit$X_sd} exists, it is used
#' automatically.
#'
#' @examples
#' \dontrun{
#' # Fit horseshoe on a design with 10 covariates and 3 treatment contrasts
#' # Design layout: [X(10), D(3), X:D(30)] = 43 columns total
#' fit <- hsreg(y, X_design, penalized = pen_flags,
#'                  n_mcmc = 1000, burnin = 500)
#'
#' # Extract CATEs at test points
#' cates <- hsreg_cates(fit, X_test = X_covariates_test,
#'                          n_x = 10, n_d = 3)
#' head(cates$cate_hat)
#' }
#'
#' @export
hsreg_cates <- function(fit, X_test, n_x, n_d,
                            X_sd = NULL, level = 0.95) {

  # -------------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------------
  if (!is.matrix(X_test) || !is.numeric(X_test)) {
    stop("'X_test' must be a numeric matrix")
  }
  if (anyNA(X_test)) {
    stop("'X_test' contains NA values. Remove or impute before extracting CATEs.")
  }

  n_test <- nrow(X_test)
  if (ncol(X_test) != n_x) {
    stop("'X_test' must have ", n_x, " columns (n_x), got ", ncol(X_test))
  }

  if (!is.numeric(n_x) || n_x < 1) stop("'n_x' must be a positive integer")
  if (!is.numeric(n_d) || n_d < 1) stop("'n_d' must be a positive integer")

  p_expected <- n_x + n_d + n_d * n_x

  # Extract beta draws
  if (inherits(fit, "hsreg")) {
    # beta_draws in an hsreg object are ALREADY rescaled to the original
    # parameterization (hsreg() divides by X_sd before storing). Do NOT
    # pick up X_sd here — that would cause double-rescaling.
    beta_draws <- fit$beta_draws   # p x n_mcmc (already on original scale)
  } else if (is.list(fit) && !is.null(fit$beta_draws)) {
    beta_draws <- fit$beta_draws
  } else if (is.matrix(fit)) {
    beta_draws <- fit
  } else {
    stop("'fit' must be a horseshoe object, a list with beta_draws, or a matrix")
  }

  p <- nrow(beta_draws)
  n_mcmc <- ncol(beta_draws)

  if (p != p_expected) {
    stop("Beta draws have ", p, " rows but expected ", p_expected,
         " (n_x + n_d + n_d*n_x = ", n_x, " + ", n_d, " + ", n_d * n_x, ")")
  }

  # -------------------------------------------------------------------------
  # Rescale beta draws if needed
  # -------------------------------------------------------------------------
  if (!is.null(X_sd)) {
    if (length(X_sd) != p) {
      stop("'X_sd' must have length ", p, ", got ", length(X_sd))
    }
    # beta_orig = beta_scaled / X_sd
    # Only rescale where X_sd != 0 and X_sd != 1
    needs_rescale <- (X_sd != 0) & (X_sd != 1)
    if (any(needs_rescale)) {
      beta_draws <- beta_draws / X_sd  # vectorized: row j divided by X_sd[j]
    }
  }

  # -------------------------------------------------------------------------
  # Compute CATEs for each contrast
  # -------------------------------------------------------------------------
  alpha <- 1 - level
  lo_prob <- alpha / 2
  hi_prob <- 1 - alpha / 2

  cate_hat <- matrix(NA_real_, n_test, n_d)
  cate_lo  <- matrix(NA_real_, n_test, n_d)
  cate_hi  <- matrix(NA_real_, n_test, n_d)

  for (w in seq_len(n_d)) {
    # gamma_w: treatment dummy coefficient (column index)
    gamma_col <- n_x + w

    # delta_w: interaction coefficients for this arm (column indices)
    delta_start <- n_x + n_d + (w - 1) * n_x + 1
    delta_end   <- delta_start + n_x - 1

    # gamma draws: n_mcmc vector (one row of beta_draws)
    gamma_draws <- beta_draws[gamma_col, ]   # n_mcmc

    # delta draws: n_x x n_mcmc matrix (subset of rows of beta_draws)
    delta_draws <- beta_draws[delta_start:delta_end, , drop = FALSE]  # n_x x n_mcmc

    # CATE draws: n_test x n_mcmc
    # tau_w^(s)(x_i) = gamma_w^(s) + x_i' delta_w^(s)
    # X_test %*% delta_draws gives n_test x n_mcmc
    # Then add gamma_w^(s) via broadcasting (each column gets + gamma_draws[s])
    tau_draws <- X_test %*% delta_draws  # n_test x n_mcmc
    tau_draws <- sweep(tau_draws, 2, gamma_draws, FUN = "+")

    # Posterior means: average across MCMC draws
    cate_hat[, w] <- rowMeans(tau_draws)

    # Credible intervals: quantiles per observation
    for (i in seq_len(n_test)) {
      qi <- hs_quantile_pair(tau_draws[i, ], lo_prob, hi_prob)
      cate_lo[i, w] <- qi["lower"]
      cate_hi[i, w] <- qi["upper"]
    }
  }

  # -------------------------------------------------------------------------
  # Return results
  # -------------------------------------------------------------------------
  list(
    cate_hat = cate_hat,
    cate_lo  = cate_lo,
    cate_hi  = cate_hi,
    n_test   = n_test,
    n_d      = n_d,
    n_draws  = n_mcmc,
    level    = level
  )
}


#' Extract CATEs from a Saved Draws File
#'
#' Loads posterior draws from a file saved by \code{hsreg(..., saving = ...)}
#' and extracts individual-level CATEs. This is useful when the original fit
#' object is not available but the draws were saved to disk.
#'
#' @param draws_file Character string. Path to the saved draws file.
#'   Supports \code{.rds} and \code{.csv} formats.
#' @param X_test Numeric matrix of dimension n_test x n_x. Test covariate
#'   values.
#' @param n_x Integer. Number of covariate columns.
#' @param n_d Integer. Number of treatment dummy columns (contrasts).
#' @param X_sd Numeric vector of length p, or \code{NULL}. Column standard
#'   deviations for rescaling. Default: NULL.
#' @param level Numeric in (0, 1). Credible interval level. Default: 0.95.
#'
#' @return Same as \code{\link{hsreg_cates}}: a list with \code{cate_hat},
#'   \code{cate_lo}, \code{cate_hi}, and metadata.
#'
#' @details
#' The draws file must have columns: \code{draw}, \code{sigma2}, \code{tau2},
#' followed by one column per beta coefficient. The beta columns are everything
#' after the first three columns.
#'
#' @seealso \code{\link{hsreg_cates}}, \code{\link{hsreg}}
#'
#' @export
hsreg_cates_from_file <- function(draws_file, X_test, n_x, n_d,
                                       X_sd = NULL, level = 0.95) {

  # -------------------------------------------------------------------------
  # Load draws
  # -------------------------------------------------------------------------
  if (!file.exists(draws_file)) {
    stop("Draws file not found: ", draws_file)
  }

  ext <- tolower(tools::file_ext(draws_file))
  if (ext == "csv") {
    draws_df <- utils::read.csv(draws_file, stringsAsFactors = FALSE)
  } else {
    # Default: .rds
    draws_df <- readRDS(draws_file)
  }

  if (!is.data.frame(draws_df)) {
    stop("Draws file must contain a data.frame, got ", class(draws_df)[1])
  }

  # The first 3 columns are draw, sigma2, tau2; the rest are beta columns
  if (ncol(draws_df) < 4) {
    stop("Draws file must have at least 4 columns (draw, sigma2, tau2, beta_1+), got ",
         ncol(draws_df))
  }

  beta_mat <- as.matrix(draws_df[, -(1:3)])  # n_mcmc x p
  # Transpose to p x n_mcmc (our convention)
  beta_draws <- t(beta_mat)

  # -------------------------------------------------------------------------
  # Call the main CATE extraction function
  # -------------------------------------------------------------------------
  hsreg_cates(
    fit = beta_draws,   # pass matrix directly
    X_test = X_test,
    n_x = n_x,
    n_d = n_d,
    X_sd = X_sd,
    level = level
  )
}
