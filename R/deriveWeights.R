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
deriveWeights <- function(norm, raw = NULL, error.model = "quantile") {
  # ----------------------------
  # Input checks
  # ----------------------------
  if (is.null(raw)) {
    stop("Raw counts must be provided for weight estimation.")
  }

  if (!is.matrix(raw)) {
    stop("raw must be a matrix (genes x samples).")
  }

  if (!identical(dim(norm), dim(raw))) {
    stop("norm and raw must have identical dimensions.")
  }

  # ----------------------------
  # Estimate technical variance
  # ----------------------------

  message("Estimating gene weights from bulk RNA-seq variance...")

  sds.tech <- runErrorModel(
    counts = raw,
    method = error.model
  )

  # ----------------------------
  # Stabilise SDs
  # ----------------------------

  sds.tech[!is.finite(sds.tech)] <- median(sds.tech, na.rm = TRUE)
  sds.tech[sds.tech < 0.05] <- 0.05
  sds.tech[sds.tech > 10]   <- 10

  # ----------------------------
  # Compute weights
  # ----------------------------

  wts <- 1 / sds.tech

  # Guard against zero / infinite weights
  wts[!is.finite(wts)] <- min(wts[is.finite(wts)], na.rm = TRUE)

  dimnames(wts) <- dimnames(norm)

  return(wts)
}
