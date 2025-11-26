#' Derive per-gene weights for BulkDecon
#'
#' Computes technical and biological variance components and returns a
#' gene-by-sample weight matrix used by \code{algorithm2()} and
#' \code{bulkdecon()}. If \code{raw} is provided, technical SDs are learned
#' from the data. Biological SDs can optionally be estimated from
#' TIL residual variances included in the package.
#'
#' @param norm Matrix of normalized expression (genes × samples).
#' @param raw Optional raw counts matrix (genes × samples).
#' @param error.model Character; currently unused except passed to
#'   \code{runErrorModel()} for future-proofing. Default: "dsp".
#' @param weight.by.TIL.resid.sd Logical; if TRUE, use the built-in
#'   \code{mean.resid.sd} vector to model biological variability.
#'
#' @return A matrix of weights of the same dimension as \code{norm}.
#'
#' @importFrom utils data
#' @export
deriveWeights <- function(norm,
                          raw = NULL,
                          error.model = "dsp",
                          weight.by.TIL.resid.sd = FALSE) {

  ## -------------------------------------------
  ## 0. Check inputs
  ## -------------------------------------------
  if (!is.matrix(norm))
    stop("`norm` must be a matrix (genes × samples).")

  G <- nrow(norm)
  N <- ncol(norm)

  ## -------------------------------------------
  ## 1. Technical SDs
  ## -------------------------------------------
  if (is.null(raw)) {

    # default small technical noise
    sds.tech <- matrix(
      0.1,
      nrow = G,
      ncol = N,
      dimnames = dimnames(norm)
    )

  } else {

    # raw provided → fit DSP error model
    if (!is.matrix(raw))
      stop("`raw` must be a matrix when provided.")

    if (!all(dim(raw) == dim(norm))) {
      warning("`raw` and `norm` have different dimensions. Attempting to align.")
      shared <- intersect(rownames(norm), rownames(raw))
      if (length(shared) == 0)
        stop("No shared genes between norm and raw.")
      raw <- raw[shared, , drop = FALSE]
      norm <- norm[shared, , drop = FALSE]
      G <- nrow(norm); N <- ncol(norm)
    }

    sds.tech <- runErrorModel(
      counts = raw,
      platform = error.model
    )

    # enforce dimnames match norm
    sds.tech <- sds.tech[rownames(norm), colnames(norm)]
  }

  ## -------------------------------------------
  ## 2. Biological SDs
  ## -------------------------------------------
  if (!weight.by.TIL.resid.sd) {

    # constant mild biological noise
    sds.bio <- matrix(
      0.1,
      nrow = G,
      ncol = N,
      dimnames = dimnames(norm)
    )

  } else {

    utils::data("mean.resid.sd", envir = environment())

    if (!exists("mean.resid.sd"))
      stop("mean.resid.sd data not found in BulkDecon package.")

    sds.bio <- matrix(
      NA,
      nrow = G,
      ncol = N,
      dimnames = dimnames(norm)
    )

    # fill in biological SDs for overlapping genes
    sharedgenes <- intersect(rownames(norm), names(mean.resid.sd))
    sds.bio[sharedgenes, ] <- mean.resid.sd[sharedgenes]

    # replace missing biological SDs with mean
    sds.bio <- replace(sds.bio, is.na(sds.bio), mean(sds.bio, na.rm = TRUE))
  }

  ## -------------------------------------------
  ## 3. Total SD and weights
  ## -------------------------------------------
  sds.tot <- sqrt(sds.tech^2 + sds.bio^2)

  # Proper inverse-variance weights
  wts <- 1 / sds.tot

  return(wts)
}
