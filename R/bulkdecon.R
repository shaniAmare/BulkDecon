#' Mixed Cell Deconvolution of Bulk Gene Expression Data
#'
#' @description
#' Performs bulk deconvolution using the BulkDecon framework, with optional
#' weighting, tumor merging, cell-type collapsing, and p-value computation.
#' This function is adapted from the SpatialDecon implementation but uses the
#' BulkDecon package structure and bundled reference matrix (safeTME).
#'
#' @details
#' The function accepts either a vector (single sample) or a matrix of bulk
#' expression data and estimates the cellular composition using a reference
#' profile matrix. Additional steps include:
#'
#' * optional gene weighting via residual stabilization
#' * merging pure-tumor samples into the reference
#' * collapsing granular cell types into user-defined merged classes
#' * computing t-statistics, standard errors, and p-values
#' * proportion scaling and optional conversion to estimated cell counts
#'
#' @param norm A p-length expression vector or a p × N matrix of log-expression.
#' @param bg Background expectation (same dimension as `norm`).
#' @param X Cell profile matrix. If `NULL`, the bundled `safeTME` reference
#'   from **BulkDecon** is used.
#' @param raw Optional raw (linear-scale) expression matrix for deriving weights.
#' @param wts Optional weights matrix for genes × samples.
#' @param resid_thresh Threshold (log2 units) for outlier residual removal.
#' @param lower_thresh Threshold for residual stabilization.
#' @param align_genes Logical; if `TRUE`, aligns the gene sets between `norm` and `X`.
#' @param is_pure_tumor Optional logical vector indicating pure tumor samples.
#' @param n_tumor_clusters Number of tumor clusters to merge into the reference.
#' @param cell_counts Optional numeric vector of nuclei counts per sample.
#' @param cellmerges Optional list mapping granular → merged cell types.
#' @param maxit Maximum number of iterations for the optimization procedure.
#'
#' @return
#' A list containing:
#' * `beta` – estimated cell-type abundance matrix
#' * `sigma` – variance estimates
#' * `p`, `t`, `se` – p-values, t-statistics, and standard errors
#' * `prop_of_all`, `prop_of_nontumor` – proportional abundances
#' * `cell.counts` – estimated cell counts (if `cell_counts` supplied)
#' * `X` – processed/merged reference matrix
#' * `beta.granular`, `sigma.granular` – pre-merge estimates (if applicable)
#' * all intermediate fields returned by the core deconvolution algorithm
#'
#' @importFrom stats pnorm
#' @importFrom utils data
#' @export
bulkdecon <- function(norm,
                      bg,
                      X = NULL,
                      raw = NULL,
                      wts = NULL,
                      resid_thresh = 3, lower_thresh = 0.5,
                      align_genes = TRUE,
                      is_pure_tumor = NULL, n_tumor_clusters = 10,
                      cell_counts = NULL,
                      cellmerges = NULL,
                      maxit = 1000) {

  #### preliminaries ---------------------------------

  # check formatting:
  if (!is.matrix(norm)) stop("norm should be a matrix")
  if ((length(X) > 0) & (!is.matrix(X))) stop("X should be a matrix")
  if ((length(raw) > 0) & (!is.matrix(raw))) stop("raw must be a matrix")
  if ((length(wts) > 0) & (!is.matrix(wts))) stop("wts must be a matrix")
  if ((length(cell_counts) > 0) & (!is.numeric(cell_counts)))
    stop("cell_counts must be numeric")

  # background handling
  if (length(bg) == 1) {
    bg <- matrix(bg, nrow(norm), ncol(norm),
                 dimnames = list(rownames(norm), colnames(norm)))
  }

  # If a matrix other than safeTME is input, rescale training matrix:
  if (length(X) > 0) {
    X <- X * 2 / quantile(X, 0.99)
  }

  # load bundled safeTME if X not supplied
  if (length(X) == 0) {
    utils::data("safeTME", package = "BulkDecon-package", envir = environment())
    X <- safeTME
  }

  sharedgenes <- intersect(rownames(norm), rownames(X))
  if (length(sharedgenes) == 0) stop("no shared gene names between norm and X")
  if (length(sharedgenes) < 100) {
    stop(paste0("Only ", length(sharedgenes),
                " genes are shared between norm and X — may not be enough."))
  }

  # calculate weights (if raw provided)
  if (length(raw) > 0) {
    weight.by.TIL.resid.sd <-
      length(intersect(colnames(X), colnames(safeTME))) > 10

    wts <- deriveWeights(
      norm,
      raw = raw,
      error.model = "dsp",
      weight.by.TIL.resid.sd = weight.by.TIL.resid.sd
    )
  }

  #### if pure tumor samples specified ------------------------------
  if (sum(is_pure_tumor) > 0) {

    X <- mergeTumorIntoX(
      norm = norm,
      bg = bg,
      pure_tumor_ids = is_pure_tumor,
      X = X[sharedgenes, ]
    )

    sharedgenes <- intersect(rownames(norm), rownames(X))
  }

  #### run decon ----------------------------------------------------
  res <- algorithm2(
    Y = norm[sharedgenes, ],
    bg = bg[sharedgenes, ],
    X = X[sharedgenes, ],
    weights = wts[sharedgenes, ],
    maxit = maxit,
    resid_thresh = resid_thresh,
    lower_thresh = lower_thresh
  )

  #### combine related cell types ----------------------------------
  if (length(cellmerges) > 0) {
    tempconv <- convertCellTypes(
      beta = res$beta,
      matching = cellmerges,
      stat = sum,
      na.rm = FALSE,
      sigma = res$sigmas
    )

    res$beta.granular <- res$beta
    res$sigma.granular <- res$sigmas

    res$beta <- tempconv$beta
    res$sigma <- tempconv$sigma
  }

  #### compute p-values --------------------------------------------
  tempbeta <- res$beta
  tempse <- tempp <- tempt <- tempbeta * NA

  for (i in seq_len(ncol(tempse))) {
    tempse[, i] <- suppressWarnings(sqrt(diag(res$sigma[, , i])))
  }

  tempt <- tempbeta / tempse
  tempp <- 2 * (1 - stats::pnorm(tempt))

  res$p <- tempp
  res$t <- tempt
  res$se <- tempse

  #### abundances: proportions + cell counts ------------------------
  res$prop_of_all <- sweep(res$beta, 2, colSums(res$beta), "/")

  nontumor <- rownames(res$beta)[!grepl("tumor", rownames(res$beta))]
  res$prop_of_nontumor <- sweep(
    res$beta[nontumor, ], 2,
    colSums(res$beta[nontumor, ]), "/"
  )

  if (length(cell_counts) > 0) {
    res$cell.counts <- convertCellScoresToCounts(
      beta = res$beta,
      nuclei.counts = cell_counts,
      omit.tumor = TRUE
    )

    if (length(res$beta.granular) > 0) {
      res$cell.counts.granular <- convertCellScoresToCounts(
        beta = res$beta.granular,
        nuclei.counts = cell_counts,
        omit.tumor = TRUE
      )
    }
  }

  res$X <- X[rownames(res$resids), ]
  return(res)
}
