#' Derive Gene Weights from Bulk RNA-seq Data
#'
#' Computes per-gene weights for deconvolution by estimating variance across
#' samples in normalized expression space. More stable genes receive higher
#' weight during regression.
#'
#' @param norm Matrix of normalized expression (genes × samples, log scale).
#' @param raw Optional raw count matrix (not used but retained for compatibility).
#' @param weight.by.TIL.resid.sd Logical; if TRUE, emphasize immune-stable genes
#'   (currently disabled for bulk RNA-seq mode).
#'
#' @return A matrix of weights (genes × samples).
#'
#' @export
deriveWeights <- function(norm,
                          raw = NULL,
                          weight.by.TIL.resid.sd = FALSE) {

  if (!is.matrix(norm))
    stop("`norm` must be a matrix (genes × samples).")

  message("Estimating gene weights from bulk RNA-seq variance...")

  ## -------------------------------------------------
  ## 1. Estimate gene-wise variance (on log scale)
  ## -------------------------------------------------
  gene_sd <- apply(norm, 1, sd, na.rm = TRUE)

  ## Guard against degenerate variance
  gene_sd[gene_sd == 0 | is.na(gene_sd)] <- median(gene_sd, na.rm = TRUE)

  ## -------------------------------------------------
  ## 2. Convert variance → weights
  ## -------------------------------------------------
  ## Lower variance → higher weight
  gene_wts <- 1 / gene_sd

  ## Normalize weights (mean = 1)
  gene_wts <- gene_wts / mean(gene_wts, na.rm = TRUE)

  ## -------------------------------------------------
  ## 3. Expand to genes × samples matrix
  ## -------------------------------------------------
  wts <- matrix(
    gene_wts,
    nrow = length(gene_wts),
    ncol = ncol(norm),
    dimnames = dimnames(norm)
  )

  return(wts)
}
