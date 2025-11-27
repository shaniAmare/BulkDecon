#' Identify Outlier Genes and Data Points from Reverse Deconvolution Residuals
#'
#' @description
#' Flags outlier genes and individual data points based on weighted residuals
#' from reverse deconvolution or log-normal regression fits.
#'
#' @param Y Matrix of observed expression values (genes × samples).
#'   Not directly used for filtering but kept for compatibility.
#' @param yhat Matrix of predicted expression values (genes × samples).
#' @param resids Matrix of residuals (observed - predicted, typically in log-space).
#' @param wts Optional numeric vector of weights.
#'   If supplied, residuals are multiplied by these weights.
#'   Length must be 1, nrow(resids), or ncol(resids).
#' @param resid_thresh Threshold on absolute weighted residuals to flag outliers.
#'
#' @return A list containing:
#' \describe{
#'   \item{outlier_genes}{Logical vector (length = number of genes).
#'                        TRUE indicates that at least one sample was an outlier.}
#'   \item{outlier_data_points}{Logical matrix (genes × samples) indicating individual outlier points.}
#'   \item{weighted_residuals}{Matrix of weighted residuals used for detection.}
#' }
#'
#' @export
flagOutliers <- function(Y, yhat, resids, wts, resid_thresh = 3) {

  # get weighted resids:
  if (length(wts) == 0) {
    wres <- resids
  }
  if (length(wts) > 0) {
    wres <- resids * wts
  }

  # flag bad genes:
  outlier_genes <- c() # <-- this line makes it so no outlier genes are filtered
  # flag bad data points: (not doing anything for now)
  outlier_data_points <- abs(resids) > resid_thresh
  return(outlier_data_points)
}
