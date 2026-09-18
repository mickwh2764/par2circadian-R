#' Fit AR(2) to every gene in an expression matrix
#'
#' @param expression A numeric matrix of genes (rows) by timepoints (columns),
#'   or anything coercible by \code{as.matrix}. Missing values are dropped per
#'   gene; genes left with fewer than 6 timepoints are skipped rather than
#'   returned with a meaningless fit.
#' @param gene_names Optional character vector of gene identifiers. Defaults to
#'   the matrix rownames, or \code{Gene_1 ... Gene_n} when there are none.
#'
#' @return A data frame of class \code{par2_batch}, one row per fitted gene,
#'   ordered by decreasing \eqn{|\lambda|}, with columns \code{gene},
#'   \code{eigenvalue}, \code{phi1}, \code{phi2}, \code{r2}, \code{root_type},
#'   \code{half_life} and \code{eigenperiod}.
#'
#' @examples
#' set.seed(1)
#' m <- t(replicate(5, as.numeric(arima.sim(list(ar = c(0.5, -0.2)), n = 24))))
#' rownames(m) <- paste0("Gene", 1:5)
#' head(fit_ar2_batch(m))
#'
#' @export
fit_ar2_batch <- function(expression, gene_names = NULL) {
  m <- as.matrix(expression)
  if (!is.numeric(m)) {
    stop("expression must be numeric", call. = FALSE)
  }
  if (is.null(gene_names)) {
    gene_names <- rownames(m)
  }
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene_", seq_len(nrow(m)))
  }
  if (length(gene_names) != nrow(m)) {
    stop("gene_names must have one entry per row of expression", call. = FALSE)
  }

  rows <- lapply(seq_len(nrow(m)), function(i) {
    values <- m[i, ]
    values <- values[!is.na(values)]
    if (length(values) < 6L) {
      return(NULL)
    }
    fit <- tryCatch(fit_ar2(values), error = function(e) NULL)
    if (is.null(fit)) {
      return(NULL)
    }
    data.frame(
      gene = gene_names[i],
      eigenvalue = fit$eigenvalue,
      phi1 = fit$phi1,
      phi2 = fit$phi2,
      r2 = fit$r2,
      root_type = fit$root_type,
      half_life = fit$half_life,
      eigenperiod = fit$eigenperiod,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(
      gene = character(0), eigenvalue = numeric(0), phi1 = numeric(0),
      phi2 = numeric(0), r2 = numeric(0), root_type = character(0),
      half_life = numeric(0), eigenperiod = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    out <- out[order(-out$eigenvalue), , drop = FALSE]
    rownames(out) <- NULL
  }
  structure(out, class = c("par2_batch", "data.frame"))
}

#' Classify a gene's dynamics from its AR(2) result
#'
#' The label depends on root type as well as magnitude: a high \eqn{|\lambda|}
#' with real roots is persistent but not oscillatory, and calling it an
#' oscillator is the most common misreading of the metric.
#'
#' @param eigenvalue Numeric \eqn{|\lambda|}.
#' @param root_type \code{"Complex"} or \code{"Real"}.
#'
#' @return One of \code{"Unstable"}, \code{"Sustained oscillator"},
#'   \code{"Damped oscillator"}, \code{"Overdamped decay"} or
#'   \code{"Rapid decay"}.
#'
#' @examples
#' classify_dynamics(0.92, "Complex")
#' classify_dynamics(0.92, "Real")
#'
#' @export
classify_dynamics <- function(eigenvalue, root_type) {
  if (eigenvalue >= 1) {
    return("Unstable")
  }
  if (identical(root_type, "Complex")) {
    if (eigenvalue >= 0.8) {
      "Sustained oscillator"
    } else if (eigenvalue >= 0.4) {
      "Damped oscillator"
    } else {
      "Rapid decay"
    }
  } else {
    if (eigenvalue >= 0.4) "Overdamped decay" else "Rapid decay"
  }
}
