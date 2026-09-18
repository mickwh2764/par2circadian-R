#' Residual-bootstrap confidence intervals for an AR(2) fit
#'
#' Resamples the fitted residuals with replacement, regenerates the series from
#' the fitted recursion, refits, and returns the percentile interval over the
#' refits. This is the standard interval for autoregressive coefficients and
#' needs only the one fitted series -- no biological replicates, which most
#' circadian designs (one animal pool per timepoint) do not provide.
#'
#' @param centred The mean-centred series the coefficients were fitted to.
#' @param phi1,phi2 Fitted AR(2) coefficients.
#' @param n_bootstrap Number of bootstrap draws.
#' @param seed Optional integer seed.
#' @param conf Confidence level as a percentage.
#'
#' @return A list with \code{eigenvalue_ci}, \code{phi1_ci}, \code{phi2_ci}
#'   (each a length-2 numeric vector) and \code{n_timepoints}.
#'
#' @section Reproducibility across implementations:
#' Point estimates from this package agree with the Python implementation to
#' within floating-point rounding, but bootstrap intervals do not: R and NumPy
#' draw from different generators, so the same seed selects different
#' resamples. Compare interval widths, not endpoints.
#'
#' @examples
#' set.seed(1)
#' x <- as.numeric(arima.sim(list(ar = c(0.6, -0.3)), n = 48))
#' fit_ar2(x, n_bootstrap = 200, seed = 42)$eigenvalue_ci
#'
#' @export
bootstrap_ar2 <- function(centred, phi1, phi2, n_bootstrap = 2000, seed = NULL,
                          conf = 95) {
  centred <- as.numeric(centred)
  n <- length(centred)
  resid <- centred[3:n] - (phi1 * centred[2:(n - 1L)] + phi2 * centred[1:(n - 2L)])
  resid <- resid - mean(resid)

  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    } else {
      NULL
    }
    on.exit(
      if (is.null(old)) {
        rm(".Random.seed", envir = globalenv())
      } else {
        assign(".Random.seed", old, envir = globalenv())
      },
      add = TRUE
    )
    set.seed(seed)
  }

  draws <- vapply(
    seq_len(n_bootstrap),
    function(i) {
      eps <- sample(resid, n - 2L, replace = TRUE)
      xb <- numeric(n)
      xb[1:2] <- centred[1:2]
      for (t in 3:n) {
        xb[t] <- phi1 * xb[t - 1L] + phi2 * xb[t - 2L] + eps[t - 2L]
      }
      refit <- tryCatch(fit_ar2(xb), error = function(e) NULL)
      if (is.null(refit)) {
        c(NA_real_, NA_real_, NA_real_)
      } else {
        c(refit$eigenvalue, refit$phi1, refit$phi2)
      }
    },
    numeric(3)
  )

  probs <- c((100 - conf) / 2, 100 - (100 - conf) / 2) / 100
  interval <- function(v) {
    v <- v[!is.na(v)]
    if (!length(v)) return(NULL)
    round(as.numeric(stats::quantile(v, probs, names = FALSE, type = 7)), 6)
  }

  list(
    eigenvalue_ci = interval(draws[1L, ]),
    phi1_ci = interval(draws[2L, ]),
    phi2_ci = interval(draws[3L, ]),
    n_timepoints = n
  )
}
