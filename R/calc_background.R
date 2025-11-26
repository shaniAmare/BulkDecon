#' Calculate Background Expression from Stable Background Genes
#'
#' This function identifies background-stable genes using `calc_bgindex()`,
#' selects high-stability genes, and uses them to estimate expected background
#' expression levels in bulk RNA-seq data. The resulting background matrix has
#' the same dimensions as `exprs_mat`.
#'
#' @param exprs_mat Numeric matrix of raw or normalized expression values.
#'        Genes in rows, samples in columns.
#' @param norm Numeric matrix of normalized expression values with the same
#'        dimensions and gene names as `exprs_mat`.
#'
#' @return A numeric matrix of background estimates with the same dimensions as
#'         `exprs_mat`.
#'
#' @details
#' Background-stable genes are defined using the `bgIdx` score produced by
#' `calc_bgindex()`. Genes ranked above the 80th percentile of the stability
#' index are considered background-like and used to compute per-sample background
#' using column means (on normalized data).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' bg <- calc_background(exprs_mat = bulk_counts, norm = normalized_counts)
#' }
calc_background <- function(exprs_mat, norm) {

  ## Check inputs
  if (!is.matrix(exprs_mat))
    stop("`exprs_mat` must be a matrix.")
  if (!is.matrix(norm))
    stop("`norm` must be a matrix.")
  if (!identical(rownames(exprs_mat), rownames(norm)))
    stop("`exprs_mat` and `norm` must have identical rownames.")

  ## Initialise output background matrix
  bg <- matrix(0,
               nrow = nrow(exprs_mat),
               ncol = ncol(exprs_mat),
               dimnames = dimnames(exprs_mat))

  ## Compute background index table
  results <- calc_bgindex(exprs_mat)

  if (!"bgIdx" %in% colnames(results)) {
    stop("`calc_bgindex()` output must contain a column named 'bgIdx'.")
  }

  ## Identify background-stable genes (top 20% bgIdx)
  cutoff <- stats::quantile(results$bgIdx, probs = 0.8, na.rm = TRUE)
  negnames <- rownames(results)[results$bgIdx > cutoff]

  ## Intersect with expressed genes
  tempnegs <- intersect(negnames, rownames(exprs_mat))

  if (length(tempnegs) == 0) {
    stop("No stable background genes identified.
              Did you supply raw or correctly normalized data?")
  }

  ## Compute per-sample background factors
  tempnegfactor <- colMeans(norm[tempnegs, , drop = FALSE])

  ## Fill background matrix
  bg <- sweep(bg, 2, tempnegfactor, "+")

  return(bg)
}
