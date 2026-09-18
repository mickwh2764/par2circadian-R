#' Persistence half-life in time units
#'
#' The time for a perturbation to decay to half its amplitude, given the
#' eigenvalue modulus and the sampling interval.
#'
#' @param eigenvalue Numeric \eqn{|\lambda|}.
#' @param sampling_interval Time between consecutive samples, in whatever units
#'   the result should be reported in. Defaults to 2 hours, the interval of the
#'   reference liver series.
#'
#' @return Half-life in the units of \code{sampling_interval}: \code{0} when
#'   \eqn{|\lambda| \le 0}, and \code{Inf} at or above 1, where the process is
#'   non-stationary and nothing decays.
#'
#' @examples
#' half_life(0.8)
#' half_life(0.8, sampling_interval = 1)
#'
#' @export
half_life <- function(eigenvalue, sampling_interval = 2) {
  if (eigenvalue <= 0) {
    return(0)
  }
  if (eigenvalue >= 1) {
    return(Inf)
  }
  sampling_interval * log(2) / (-log(eigenvalue))
}

#' Intrinsic eigenperiod in time units
#'
#' The oscillation period implied by a complex-root AR(2) fit. This is the
#' period the dynamics prefer, which need not equal the period of any imposed
#' driver.
#'
#' @param phi1,phi2 AR(2) coefficients.
#' @param sampling_interval Time between consecutive samples.
#'
#' @return The period in the units of \code{sampling_interval}, or \code{NA} for
#'   real roots, where there is no oscillation to have a period.
#'
#' @examples
#' eigenperiod(0.6, -0.3)
#' eigenperiod(1.2, -0.2)
#'
#' @export
eigenperiod <- function(phi1, phi2, sampling_interval = 2) {
  disc <- phi1^2 + 4 * phi2
  if (disc >= 0) {
    return(NA_real_)
  }
  omega <- atan2(sqrt(-disc), phi1)
  if (omega <= 0) {
    return(NA_real_)
  }
  sampling_interval * 2 * pi / omega
}
