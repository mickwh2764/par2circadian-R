#' Fit an AR(2) model to one expression time series
#'
#' Fits \eqn{x_t = \phi_1 x_{t-1} + \phi_2 x_{t-2} + \epsilon_t} to a
#' mean-centred series by ordinary least squares and returns the eigenvalue
#' modulus \eqn{|\lambda|} of the characteristic equation
#' \eqn{r^2 - \phi_1 r - \phi_2 = 0}, together with the derived persistence
#' measures.
#'
#' \eqn{|\lambda|} measures how strongly recent history constrains future
#' values. It is not a test of rhythmicity: a non-oscillatory red-noise series
#' can score high, correctly, and a clean fixed-period cosine need not. Use it
#' alongside a rhythm-detection method, not instead of one.
#'
#' Roots are complex when \eqn{\phi_1^2 + 4\phi_2 < 0}, which indicates
#' oscillatory dynamics; real roots indicate overdamped, non-oscillatory decay.
#'
#' @param x Numeric vector of expression values, at least 6 timepoints.
#' @param n_bootstrap Number of residual-bootstrap draws. When greater than
#'   zero, percentile confidence intervals are added; 2000 is ample.
#' @param seed Optional integer seed for reproducible bootstrap intervals.
#' @param conf Confidence level for the bootstrap intervals, as a percentage.
#'
#' @return An object of class \code{par2_fit}: a list with elements
#'   \code{eigenvalue} (\eqn{|\lambda|}), \code{phi1}, \code{phi2}, \code{r2},
#'   \code{root_type} (\code{"Complex"} or \code{"Real"}), \code{half_life} in
#'   samples, and \code{eigenperiod} in samples (\code{NA} for real roots).
#'   With \code{n_bootstrap > 0} it also carries \code{eigenvalue_ci},
#'   \code{phi1_ci}, \code{phi2_ci} and \code{n_timepoints}.
#'
#' @section Precision:
#' A point estimate from a short series is imprecise: at 24 evenly sampled
#' timepoints the 95% interval is typically about 0.4 wide, and the estimator
#' is biased upward below roughly 24 points. Report the interval alongside any
#' per-gene value and do not interpret differences smaller than it.
#'
#' @examples
#' set.seed(1)
#' x <- as.numeric(arima.sim(list(ar = c(0.6, -0.3)), n = 48))
#' fit <- fit_ar2(x)
#' fit$eigenvalue
#' fit$root_type
#'
#' @seealso \code{\link{fit_ar2_batch}} for a whole expression matrix,
#'   \code{\link{classify_dynamics}}, \code{\link{half_life}}.
#' @export
fit_ar2 <- function(x, n_bootstrap = 0, seed = NULL, conf = 95) {
  x <- as.numeric(x)
  if (length(x) < 6L) {
    stop(sprintf("Need >= 6 timepoints, got %d", length(x)), call. = FALSE)
  }
  if (anyNA(x)) {
    stop("x must not contain missing values; drop or impute them first", call. = FALSE)
  }

  centred <- x - mean(x)
  n <- length(centred)

  y <- centred[3:n]
  design <- cbind(centred[2:(n - 1L)], centred[1:(n - 2L)])

  phi <- ar2_coefficients(design, y)
  phi1 <- phi[1L]
  phi2 <- phi[2L]

  fitted <- as.numeric(design %*% phi)
  ss_res <- sum((y - fitted)^2)
  ss_tot <- sum((y - mean(y))^2)
  r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else 0

  roots <- eigen_from_coefficients(phi1, phi2)

  result <- list(
    eigenvalue = round(roots$eigenvalue, 6),
    phi1 = round(phi1, 6),
    phi2 = round(phi2, 6),
    r2 = round(r2, 6),
    root_type = roots$root_type,
    half_life = persistence_half_life(roots$eigenvalue),
    eigenperiod = if (is.na(roots$eigenperiod)) NA_real_ else round(roots$eigenperiod, 4)
  )

  if (n_bootstrap > 0) {
    result <- c(result, bootstrap_ar2(centred, phi1, phi2,
      n_bootstrap = n_bootstrap, seed = seed, conf = conf
    ))
  }

  structure(result, class = "par2_fit")
}

#' @keywords internal
#' @noRd
ar2_coefficients <- function(design, y) {
  xtx <- crossprod(design)
  xty <- crossprod(design, y)
  phi <- tryCatch(solve(xtx, xty),
    error = function(e) qr.coef(qr(design), y)
  )
  as.numeric(phi)
}

#' @keywords internal
#' @noRd
eigen_from_coefficients <- function(phi1, phi2) {
  disc <- phi1^2 + 4 * phi2
  if (disc < 0) {
    eigenvalue <- sqrt(-phi2)
    omega <- atan2(sqrt(-disc), phi1)
    eigenperiod <- if (omega > 0) 2 * pi / omega else NA_real_
    root_type <- "Complex"
  } else {
    sqrt_disc <- sqrt(disc)
    eigenvalue <- max(abs((phi1 + sqrt_disc) / 2), abs((phi1 - sqrt_disc) / 2))
    eigenperiod <- NA_real_
    root_type <- "Real"
  }
  list(
    eigenvalue = min(eigenvalue, 2),
    eigenperiod = eigenperiod,
    root_type = root_type
  )
}

#' @keywords internal
#' @noRd
persistence_half_life <- function(eigenvalue) {
  if (eigenvalue > 0 && eigenvalue < 1) {
    round(log(2) / (-log(eigenvalue)), 4)
  } else {
    NA_real_
  }
}

#' @param ... Ignored.
#' @rdname fit_ar2
#' @export
print.par2_fit <- function(x, ...) {
  cat(sprintf("AR(2) fit: |lambda| = %.6f (%s roots)\n", x$eigenvalue, tolower(x$root_type)))
  if (!is.null(x$eigenvalue_ci)) {
    cat(sprintf(
      "  bootstrap CI  [%.6f, %.6f] from %d timepoints\n",
      x$eigenvalue_ci[1L], x$eigenvalue_ci[2L], x$n_timepoints
    ))
  }
  cat(sprintf("  phi1 = %.6f, phi2 = %.6f, R^2 = %.6f\n", x$phi1, x$phi2, x$r2))
  if (!is.na(x$half_life)) {
    cat(sprintf("  half-life = %.4f samples\n", x$half_life))
  }
  if (!is.na(x$eigenperiod)) {
    cat(sprintf("  eigenperiod = %.4f samples\n", x$eigenperiod))
  }
  invisible(x)
}
