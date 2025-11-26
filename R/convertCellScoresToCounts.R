#' Convert Cell Abundance Scores to Estimated Counts
#'
#' This function rescales cell-type abundance scores (e.g., from
#' bulk RNA-seq deconvolution) into interpretable units:
#' *cells per 100 nuclei* and optionally *absolute cell counts*
#' if nuclei counts per sample are provided.
#'
#' @param beta A matrix or data.frame of cell-type abundance scores
#'   (rows = cell types, columns = samples).
#' @param nuclei.counts Optional numeric vector giving the total nuclei
#'   counts per sample. Must be the same length as the number of columns
#'   in \code{beta}.
#' @param omit.tumor Logical; if TRUE, removes rows whose rownames
#'   contain the string `"tumor"`.
#'
#' @details
#' The function rescales abundance scores so that the sample with the highest
#' total abundance sums to 100 "pseudo-cells". Other samples are scaled
#' proportionally.
#'
#' If \code{nuclei.counts} is supplied, the rescaled values are converted
#' to estimated absolute cell counts.
#'
#' @return A list containing:
#' \describe{
#'   \item{cells.per.100}{A matrix of relative abundances (cells per 100 nuclei).}
#'   \item{cell.counts}{(Optional) A matrix of absolute estimated cell counts.}
#' }
#'
#' @export
convertCellScoresToCounts <- function(beta,
                                      nuclei.counts = NULL,
                                      omit.tumor = FALSE) {

  # ---- Input checks -------------------------------------------------------

  if (is.null(beta)) {
    stop("`beta` must be supplied.")
  }
  beta <- as.matrix(beta)

  if (!is.null(nuclei.counts)) {
    if (length(nuclei.counts) != ncol(beta)) {
      stop("Length of `nuclei.counts` must match number of samples (ncol(beta)).")
    }
    if (any(!is.finite(nuclei.counts)) || any(nuclei.counts <= 0)) {
      stop("`nuclei.counts` must contain positive numeric values.")
    }
  }

  # ---- Optional removal of tumor rows -------------------------------------

  if (omit.tumor) {
    beta <- beta[!grepl("tumor", rownames(beta), ignore.case = TRUE), , drop = FALSE]
  }

  # ---- Rescaling -----------------------------------------------------------

  max.total <- max(colSums(beta), na.rm = TRUE)

  if (!is.finite(max.total) || max.total == 0) {
    stop("Sum of abundance scores is zero or invalid; cannot rescale.")
  }

  cells.per.100 <- beta / max.total * 100

  # ---- Construct output ----------------------------------------------------

  out <- list(cells.per.100 = cells.per.100)

  if (!is.null(nuclei.counts)) {
    cell.counts <- sweep(cells.per.100, 2, nuclei.counts, "*") / 100
    out$cell.counts <- cell.counts
  }

  return(out)
}
