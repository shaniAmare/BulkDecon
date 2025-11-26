#' Calculate Background Expression from Stable Background Genes
#'
#' This function identifies background-stable genes using `calc_bgindex()`,
#' selects high-stability genes, and uses them to estimate expected background
#' expression levels in bulk RNA-seq data. The resulting background matrix has
#' the same dimensions as `raw`.
#'
#' @param raw Numeric matrix of *raw* bulk expression values.
#'        Genes in rows, samples in columns.
#' @param norm Numeric matrix of *normalized* expression values
#'        with the same dimensions and gene names as `raw`.
#'
#' @return A numeric matrix of background expression estimates with the same
#'         dimensions and row/column names as `raw`.
#'
#' @details
#' Background-stable genes are defined using the `bgIdx` score produced by
#' `calc_bgindex()`. Genes ranked above the 80th percentile of the stability
#' index are considered background-like and used to compute per-sample background
#' levels using column means of the normalized expression matrix (`norm`).
#'
#' These estimated background levels are then broadcast across all rows (genes)
#' to form a background matrix of the same shape as `raw`.
#'
#' @seealso
#' `calc_bgindex()` for computing background gene stability indices.
#'
#' @examples
#' \dontrun{
#' bg <- calc_background(raw = bulk_counts,
#'                       norm = normalized_counts)
#' }
#'
#' @export
calc_background <- function(raw, norm) {

  ## Check inputs
  if (!is.matrix(raw))
    stop("`raw` must be a matrix.")
  if (!is.matrix(norm))
    stop("`norm` must be a matrix.")
  if (!identical(rownames(raw), rownames(norm)))
    stop("`raw` and `norm` must have identical rownames.")

  ## Initialise output background matrix
  bg <- matrix(
    0,
    nrow = nrow(raw),
    ncol = ncol(raw),
    dimnames = dimnames(raw)
  )

  ## Compute background index table
  results <- calc_bgindex(raw)

  if (!"bgIdx" %in% colnames(results)) {
    stop("`calc_bgindex()` output must contain a column named 'bgIdx'.")
  }

  ## Identify background-stable genes (top 20% bgIdx)
  cutoff <- stats::quantile(results$bgIdx, probs = 0.8, na.rm = TRUE)
  negnames <- rownames(results)[results$bgIdx > cutoff]

  ## Intersect with expressed genes
  tempnegs <- intersect(negnames, rownames(raw))

  if (length(tempnegs) == 0) {
    stop("No stable background genes identified.
          Did you supply raw or normalized data correctly?")
  }

  ## Compute per-sample background factors
  tempnegfactor <- colMeans(norm[tempnegs, , drop = FALSE])

  ## Fill background matrix
  bg <- sweep(bg, 2, tempnegfactor, "+")

  return(bg)
}
