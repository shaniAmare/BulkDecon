#' Run the BulkDecon workflow on bulk expression data
#'
#' This function runs the full BulkDecon pipeline using:
#' \itemize{
#'   \item a normalized bulk expression matrix (`norm_elt`)
#'   \item a raw (unnormalized) bulk expression matrix (`raw_elt`)
#'   \item a single-cell or custom signature matrix (`X`)
#' }
#'
#' The function automatically:
#' \enumerate{
#'   \item converts inputs to matrices
#'   \item estimates background using \code{calc_background()}
#'   \item runs \code{bulkdecon()} with the appropriate inputs
#' }
#'
#' @param X Signature matrix (genes × cell types).
#' @param norm_elt Normalized bulk expression matrix.
#' @param raw_elt Raw (unnormalized) bulk expression matrix.
#' @param wts Optional weights matrix (same dimensions as norm_elt).
#' @param resid_thresh Residual threshold for outlier filtering (log2 scale).
#' @param lower_thresh Lower bound for stabilized log2 residuals.
#' @param align_genes Logical; align genes by shared rownames (default TRUE).
#' @param is_pure_tumor Optional logical vector indicating tumor-pure samples.
#' @param n_tumor_clusters Integer; number of tumor clusters for pure-tumor mode.
#' @param cell_counts Optional vector of known cell counts per sample.
#' @param cellmerges Optional list describing merged cell type groups.
#' @param maxit Maximum number of iterations for deconvolution.
#'
#' @return A list returned by \code{bulkdecon()}, containing:
#' \itemize{
#'   \item \code{beta}: estimated cell proportions
#'   \item \code{sigmas}: covariance matrices
#'   \item \code{p}: p-values
#'   \item \code{t}: t-statistics
#'   \item \code{se}: standard errors
#'   \item \code{yhat}: fitted values
#'   \item \code{resids}: residual matrix
#' }
#'
#' @export
runbulkdecon <- function(
    X = NULL,
    norm_elt = NULL,
    raw_elt = NULL,
    wts = NULL,
    resid_thresh = 3,
    lower_thresh = 0.5,
    align_genes = TRUE,
    is_pure_tumor = NULL,
    n_tumor_clusters = 10,
    cell_counts = NULL,
    cellmerges = NULL,
    maxit = 1000
) {

  ## ---------------------------------------------------------------
  ## 1. Validate inputs
  ## ---------------------------------------------------------------
  if (is.null(norm_elt)) {
    stop("`norm_elt` must be provided.")
  }
  if (is.null(raw_elt)) {
    stop("`raw_elt` must be provided.")
  }

  ## Convert inputs to numeric matrices
  norm <- tryCatch(as.matrix(norm_elt), error = function(e)
    stop("`norm_elt` cannot be converted to a numeric matrix.")
  )

  raw <- tryCatch(as.matrix(raw_elt), error = function(e)
    stop("`raw_elt` cannot be converted to a numeric matrix.")
  )

  ## Validate rownames
  if (is.null(rownames(norm)) || is.null(rownames(raw))) {
    stop("Both `norm_elt` and `raw_elt` must have gene rownames.")
  }

  ## ---------------------------------------------------------------
  ## 2. Estimate background – using your BulkDecon calc_background()
  ## ---------------------------------------------------------------
  bg <- calc_background(raw = raw, norm = norm)

  ## ---------------------------------------------------------------
  ## 3. Run BulkDecon core algorithm
  ## ---------------------------------------------------------------
  result <- bulkdecon(
    norm = norm,
    bg = bg,
    X = X,
    raw = raw,
    wts = wts,
    resid_thresh = resid_thresh,
    lower_thresh = lower_thresh,
    align_genes = align_genes,
    is_pure_tumor = is_pure_tumor,
    n_tumor_clusters = n_tumor_clusters,
    cell_counts = cell_counts,
    cellmerges = cellmerges,
    maxit = maxit
  )

  return(result)
}

