#' Merge Tumor Profiles Into a Reference Matrix
#'
#' @description
#' This function takes a normalized expression matrix (`norm`) and its
#' corresponding background matrix (`bg`), isolates pure tumor AOIs,
#' background-subtracts them, rescales tumor profiles to match the dynamic
#' range of the reference matrix `X`, and appends these tumor profiles to `X`.
#'
#' The function is adapted from the Nanostring SpatialDecon workflow, but
#' streamlined for bulk RNA-seq deconvolution use cases.
#'
#' @param norm A numeric matrix of normalized expression values. Rows = genes,
#'   columns = AOIs or samples.
#' @param bg A numeric matrix of background estimates, same dimensions as `norm`.
#' @param pure_tumor_ids A character vector of column names (or indices)
#'   indicating tumor-only samples to use in constructing tumor profiles.
#' @param X A reference profile matrix (e.g., cell type signatures) with genes
#'   as rows and reference profiles as columns.
#'
#' @details
#' Steps performed:
#' \enumerate{
#'   \item Replace zeros with the minimum nonzero value.
#'   \item Subset `norm` and `bg` to pure tumor AOIs.
#'   \item Background-subtract tumor expression.
#'   \item Align genes shared between tumor profiles and reference matrix `X`.
#'   \item Rescale tumor profiles to have similar 90th-percentile expression as `X`.
#'   \item Append tumor profiles to the reference matrix.
#' }
#'
#' @return A matrix combining the original reference profiles and new tumor
#'   profiles. Gene order is preserved among shared genes.
#'
#' @examples
#' \dontrun{
#' merged <- mergeTumorIntoX(norm = tumor_norm,
#'                           bg = tumor_bg,
#'                           pure_tumor_ids = c("Tumor1","Tumor2"),
#'                           X = reference_profiles)
#' }
#'
#' @export
mergeTumorIntoX <- function(norm, bg, pure_tumor_ids, X) {

  # round up 0 values in norm:
  min.nonzero <- min(norm[norm > 0], na.rm = TRUE)
  norm <- pmax(norm, min.nonzero)

  # subset data to only the pure tumor IDs:
  norm <- norm[, pure_tumor_ids, drop = FALSE]
  bg <- bg[, pure_tumor_ids, drop = FALSE]

  # background subtraction:
  norm <- pmax(norm - bg, min(norm) / 20)

  # align genes in norm with reference matrix X:
  sharedgenes <- intersect(rownames(norm), rownames(X))
  tumorX <- norm[sharedgenes, , drop = FALSE]
  X <- X[sharedgenes, , drop = FALSE]

  # rescale tumor profiles by 90th percentile:
  meanq90 <- max(mean(apply(X, 2, stats::quantile, 0.9)), 1e-3)
  tumorq90s <- apply(tumorX, 2, stats::quantile, 0.9)
  tumorX <- sweep(tumorX, 2, tumorq90s, "/") * meanq90

  # merge tumor profiles into reference X:
  out <- cbind(X, tumorX)
  return(out)
}
