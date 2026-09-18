#' Core clock and known circadian target gene symbols
#'
#' Mouse and human symbols for the transcriptional-translational feedback loop
#' (\code{CORE_CLOCK_GENES}) and for well-established clock-controlled targets
#' (\code{KNOWN_TARGET_GENES}), used to stratify a batch fit into layers.
#' Matching is by exact symbol, so supply your own set for other species or
#' nomenclatures.
#'
#' @format Character vectors.
#' @name gene_sets
#' @export
CORE_CLOCK_GENES <- c(
  "Per1", "Per2", "Per3", "Cry1", "Cry2", "Clock", "Arntl", "Bmal1",
  "Nr1d1", "Nr1d2", "Dbp", "Tef", "Npas2", "Rorc", "Rora",
  "ARNTL", "PER1", "PER2", "PER3", "CRY1", "CRY2", "CLOCK",
  "NR1D1", "NR1D2", "DBP", "TEF", "NPAS2", "RORC", "RORA"
)

#' @rdname gene_sets
#' @export
KNOWN_TARGET_GENES <- c(
  "Myc", "Ccnd1", "Wee1", "Chek2", "Tp53", "Cdkn1a", "Bcl2", "Bax",
  "Ccne1", "Ccne2", "Mcm6", "Mki67", "Lgr5", "Axin2",
  "MYC", "CCND1", "WEE1", "CHEK2", "TP53", "CDKN1A", "BCL2", "BAX",
  "CCNE1", "CCNE2", "MCM6", "MKI67", "LGR5", "AXIN2"
)

#' Assign a gene to a hierarchy layer
#'
#' @param gene_name Gene symbol.
#'
#' @return \code{"Clock"}, \code{"Target"} or \code{"Background"}.
#'
#' @examples
#' classify_gene_layer("Per2")
#' classify_gene_layer("Actb")
#'
#' @export
classify_gene_layer <- function(gene_name) {
  if (gene_name %in% CORE_CLOCK_GENES) {
    "Clock"
  } else if (gene_name %in% KNOWN_TARGET_GENES) {
    "Target"
  } else {
    "Background"
  }
}

#' Recover the three-layer persistence hierarchy from a batch fit
#'
#' Stratifies fitted genes into clock, target and background layers and reports
#' their median \eqn{|\lambda|}. The clock-to-target difference is the gearbox
#' gap: in intact tissue clock genes are the most persistent layer, and the gap
#' collapses when the feedback loop is disrupted.
#'
#' @param results A \code{par2_batch} data frame from
#'   \code{\link{fit_ar2_batch}}, or any data frame with \code{gene} and
#'   \code{eigenvalue} columns.
#' @param clock_genes Optional character vector replacing
#'   \code{CORE_CLOCK_GENES}, for other species or custom definitions.
#' @param target_genes Optional character vector replacing
#'   \code{KNOWN_TARGET_GENES}.
#'
#' @return A list with the three layer medians, \code{gearbox_gap},
#'   \code{hierarchy_preserved}, the layer sizes, the clock and target genes
#'   with their eigenvalues, and \code{health_grade}.
#'
#' @section Interpreting the grade:
#' The A-F grade is a convenience banding of the gap (A at 0.15 and above,
#' then 0.10, 0.05, 0.02, F below). The thresholds are descriptive, taken from
#' the ranges seen in intact mouse liver, and are not a validated diagnostic.
#'
#' @examples
#' set.seed(1)
#' m <- t(replicate(6, as.numeric(arima.sim(list(ar = c(0.5, -0.2)), n = 24))))
#' rownames(m) <- c("Per2", "Arntl", "Nr1d1", "Myc", "Wee1", "Actb")
#' discover_hierarchy(fit_ar2_batch(m))$gearbox_gap
#'
#' @export
discover_hierarchy <- function(results, clock_genes = NULL, target_genes = NULL) {
  if (!all(c("gene", "eigenvalue") %in% names(results))) {
    stop("results needs 'gene' and 'eigenvalue' columns", call. = FALSE)
  }
  if (is.null(clock_genes)) clock_genes <- CORE_CLOCK_GENES
  if (is.null(target_genes)) target_genes <- KNOWN_TARGET_GENES

  layer <- ifelse(results$gene %in% clock_genes, "Clock",
    ifelse(results$gene %in% target_genes, "Target", "Background")
  )
  clock <- results[layer == "Clock", c("gene", "eigenvalue")]
  target <- results[layer == "Target", c("gene", "eigenvalue")]
  background <- results$eigenvalue[layer == "Background"]

  layer_median <- function(v) if (length(v)) stats::median(v) else 0
  clock_med <- layer_median(clock$eigenvalue)
  target_med <- layer_median(target$eigenvalue)
  background_med <- layer_median(background)
  gap <- clock_med - target_med

  grade <- if (gap >= 0.15) {
    "A"
  } else if (gap >= 0.10) {
    "B"
  } else if (gap >= 0.05) {
    "C"
  } else if (gap >= 0.02) {
    "D"
  } else {
    "F"
  }

  order_by_eigenvalue <- function(df) {
    df <- df[order(-df$eigenvalue), , drop = FALSE]
    rownames(df) <- NULL
    df
  }

  list(
    clock_median = round(clock_med, 4),
    target_median = round(target_med, 4),
    background_median = round(background_med, 4),
    gearbox_gap = round(gap, 4),
    hierarchy_preserved = clock_med > target_med && target_med > background_med,
    n_clock = nrow(clock),
    n_target = nrow(target),
    n_background = length(background),
    clock_genes = order_by_eigenvalue(clock),
    target_genes = order_by_eigenvalue(target),
    health_grade = grade
  )
}

#' Gearbox gap between two sets of eigenvalues
#'
#' @param clock_eigenvalues,target_eigenvalues Numeric vectors of
#'   \eqn{|\lambda|}.
#'
#' @return \code{median(clock) - median(target)}, or \code{0} if either set is
#'   empty.
#'
#' @examples
#' gearbox_gap(c(0.72, 0.70), c(0.61, 0.59))
#'
#' @export
gearbox_gap <- function(clock_eigenvalues, target_eigenvalues) {
  if (!length(clock_eigenvalues) || !length(target_eigenvalues)) {
    return(0)
  }
  stats::median(clock_eigenvalues) - stats::median(target_eigenvalues)
}
