#' Derive Gene Weights from Bulk RNA-seq Data
#'
#' Computes per-gene weights for deconvolution by estimating variance across
#' samples in normalized expression space. More stable genes receive higher
#' weight during regression.
#'
#' @param norm Normalized expression matrix (genes x samples)
#' @param raw Raw count matrix (genes x samples): optional
#' @param error.model Variance model: "quantile" or "loglinear"
#'
#' @return A matrix of weights (genes × samples).
#'
#' @export
deriveWeights <- function(norm, raw = NULL, error.model = "bulk") {

  # get tech SDs if raw data provided:
  if (length(raw) == 0) {
    sds.tech <- matrix(0.1, nrow(raw), ncol(raw), dimnames = dimnames(raw))
  }
  if (length(raw) > 0) {
    sds.tech <- runErrorModel(
      counts = raw,
      platform = error.model
    )
  }

  # if the mean.resid.sd vector (which defines genes' biological SD) is in
  # the environment, get biological noise:
  sds.bio <- matrix(0.1, nrow(raw), ncol(raw), dimnames = dimnames(raw))

  # define total SD, and invert to get weights
  sds.tot <- sqrt(sds.tech^2 + sds.bio^2)
  wts <- 1 / sds.tech
  return(wts)
}
