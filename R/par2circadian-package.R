#' par2circadian: AR(2) eigenvalue persistence for circadian time series
#'
#' Fits a second-order autoregressive model per gene and reports the eigenvalue
#' modulus \eqn{|\lambda|} of its characteristic equation as a scalar measure of
#' temporal persistence.
#'
#' Start with \code{\link{fit_ar2}} for one series or \code{\link{fit_ar2_batch}}
#' for a matrix, then \code{\link{discover_hierarchy}} to stratify the result by
#' gene layer. \code{vignette("par2circadian")} works through a mouse liver
#' series end to end.
#'
#' \strong{What the metric is not.} \eqn{|\lambda|} does not test for
#' rhythmicity. Red noise with no rhythm at all scores high, correctly, and a
#' clean fixed-period cosine need not. On a benchmark of noise-driven
#' oscillators sampled at 24 timepoints it detected 92\% against 40-62\% for
#' JTK_CYCLE, ARSER and RAIN; on a clean 24 h cosine those methods reached 100\%
#' and \eqn{|\lambda|} reached 12\%. It answers a different question, so run it
#' alongside a rhythm-detection method rather than in place of one.
#'
#' Results agree with the reference Python implementation par2-circadian: the
#' test suite refits a 300-gene mouse liver dataset and every gene matches
#' exactly at the precision both implementations report.
#'
#' @section Citation:
#' Cite the method preprint (\doi{10.21203/rs.3.rs-9283100/v1}) alongside this
#' package; \code{citation("par2circadian")} gives both.
#'
#' @keywords internal
"_PACKAGE"
