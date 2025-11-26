#' Create a Cell-Type Profile Matrix from Single-Cell Expression Data
#'
#' This function generates a "reference profile matrix" (mean expression per
#' cell type) from single-cell count data, analogous to reference matrices used
#' for deconvolution methods. It supports filtering low-quality cells,
#' normalisation, optional gene list restriction, and saving the resulting
#' matrix to file.
#'
#' @param mtx A gene-by-cell matrix of raw counts (dense or sparse).
#' @param cellAnnots A data frame containing cell annotations.
#' @param cellTypeCol Name of the column in \code{cellAnnots} containing
#'   cell-type labels.
#' @param cellNameCol Name of the column in \code{cellAnnots} containing
#'   cell barcodes / IDs.
#' @param matrixName Base name for output file (default: "Custom").
#' @param outDir Directory for writing output (default: current directory).
#'   If \code{NULL}, no file is written.
#' @param geneList Optional character vector of genes to retain.
#' @param normalize Logical; if TRUE, performs median normalisation.
#' @param scalingFactor Numeric scale multiplier applied to final matrix.
#' @param minCellNum Minimum number of viable cells per cell type.
#' @param minGenes Minimum genes expressed per cell (non-zero counts).
#' @param discardCellTypes Logical; if TRUE, removes unknown / doublet /
#'   filtered cell types.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Validates input, ensuring cell IDs match matrix column names.
#'   \item Optionally removes ambiguous cell types.
#'   \item Optionally normalises counts by median library size.
#'   \item For each cell type:
#'     \item Filters low-quality cells.
#'     \item Computes mean expression (if enough cells pass filters).
#'   \item Filters genes expressed in at least one cell type.
#'   \item Applies scaling and optional gene restriction.
#'   \item Writes a CSV file if \code{outDir} is provided.
#' }
#'
#' @return A numeric gene-by-celltype matrix of averaged expression.
#' @export
create_profile_matrix <- function(mtx, cellAnnots, cellTypeCol, cellNameCol,
                                  matrixName = "Custom", outDir = "./",
                                  geneList = NULL, normalize = FALSE,
                                  scalingFactor = 5, minCellNum = 15,
                                  minGenes = 100, discardCellTypes = FALSE) {

  # ----------------------------- Input checks ------------------------------
  if (is.null(mtx))
    stop("A count matrix (`mtx`) must be supplied.")

  if (!is.null(outDir) && !dir.exists(outDir))
    stop("Output directory does not exist: ", outDir)

  if (is.null(cellAnnots) || !is.data.frame(cellAnnots))
    stop("`cellAnnots` must be a data frame.")

  if (is.null(cellTypeCol) || !cellTypeCol %in% colnames(cellAnnots))
    stop("`cellTypeCol` must be a valid column in `cellAnnots`.")

  if (is.null(cellNameCol) || !cellNameCol %in% colnames(cellAnnots))
    stop("`cellNameCol` must be a valid column in `cellAnnots`.")

  if (!is.logical(normalize))
    warning("`normalize` should be TRUE/FALSE; continuing with default FALSE.")

  if (!is.logical(discardCellTypes))
    warning("`discardCellTypes` should be TRUE/FALSE; continuing with default TRUE.")

  if (!is.numeric(scalingFactor)) {
    warning("`scalingFactor` must be numeric. Using 5.")
    scalingFactor <- 5
  }

  if (!is.numeric(minCellNum)) {
    warning("`minCellNum` must be numeric. Using 15.")
    minCellNum <- 15
  }

  if (!is.numeric(minGenes)) {
    warning("`minGenes` must be numeric. Using 100.")
    minGenes <- 100
  }

  # Ensure sparse
  mtx <- Matrix::Matrix(as.matrix(mtx), sparse = TRUE)

  # ----------------------- Extract cell-type annotation --------------------
  cellTypes <- cellAnnots[[cellTypeCol]]
  names(cellTypes) <- cellAnnots[[cellNameCol]]

  if (is.null(cellTypes))
    stop("Could not extract cell types from `cellAnnots`.")

  # orientation check: transpose if needed
  if (!any(names(cellTypes) %in% colnames(mtx)) &&
      any(names(cellTypes) %in% rownames(mtx))) {
    message("Transposing matrix to match cell orientation.")
    mtx <- t(mtx)
  }

  if (!any(names(cellTypes) %in% colnames(mtx)))
    stop("Cell IDs do not match any column names in `mtx`.")

  if (!all(names(cellTypes) %in% colnames(mtx))) {
    missing <- sum(!names(cellTypes) %in% colnames(mtx))
    warning(missing, " cells in cellAnnots not present in the matrix.")
  }

  # ------------------------ Remove ambiguous cell types --------------------
  if (discardCellTypes) {
    bad <- is.na(cellTypes) |
      tolower(cellTypes) %in% c("unknown", "unspecified", "not available")

    # filter common unwanted labels
    bad2 <- grep("doublet|dividing|low q|filtered|mitotic",
                 tolower(cellTypes))

    rmCells <- unique(c(which(bad), bad2))
    if (length(rmCells) > 0)
      cellTypes <- cellTypes[-rmCells]
  }

  # ---------------------------- Normalisation ------------------------------
  if (normalize) {
    message("Normalizing matrix by median library size.")
    med <- median(Matrix::colSums(mtx))
    mtx <- sweep(mtx, 2, Matrix::colSums(mtx), "/")
    mtx <- mtx * med
  }

  # ---------------------------- Build Atlas --------------------------------
  CTs <- unique(cellTypes)
  atlas <- NULL

  message("Creating profile matrix:")

  for (ct in CTs) {

    message(which(CTs == ct), "/", length(CTs), ": ", ct)

    cells <- names(cellTypes)[cellTypes == ct]
    cells <- cells[cells %in% colnames(mtx)]

    if (length(cells) < minCellNum) {
      warning("Skipping '", ct,
              "' (< minCellNum cells). Adjust thresholds if needed.")
      next
    }

    # Filter cells with low total detected genes
    detected <- Matrix::colSums(mtx[, cells] > 0)
    cells <- cells[detected > minGenes]

    if (length(cells) < minCellNum) {
      warning("Skipping '", ct,
              "' after filtering for minGenes.")
      next
    }

    # Compute mean expression for that cell type
    if (length(cells) > 1) {
      avg <- Matrix::rowMeans(mtx[, cells, drop = FALSE])
    } else {
      avg <- mtx[, cells, drop = FALSE]
    }

    # Append
    atlas <- cbind(atlas, avg)
    colnames(atlas)[ncol(atlas)] <- ct
  }

  if (is.null(atlas))
    stop("No cell types passed thresholds — profile matrix empty.")

  # ----------------------------- Gene filtering ----------------------------
  keepGenes <- Matrix::rowSums(atlas > 0) >= 1
  atlas <- atlas[keepGenes, ]

  # ------------------------------ Scaling ----------------------------------
  atlas <- atlas * scalingFactor

  # ------------------------------- Gene list --------------------------------
  if (!is.null(geneList)) {
    keep <- intersect(rownames(atlas), geneList)
    if (length(keep) == 0) {
      warning("No geneList genes found in atlas — no filtering performed.")
    } else {
      atlas <- atlas[keep, , drop = FALSE]
    }
  }

  # ------------------------------ Output -----------------------------------
  colnames(atlas) <- gsub(",", "-", colnames(atlas))

  if (!is.null(outDir)) {
    outfile <- file.path(outDir, paste0(matrixName, "_profileMatrix.csv"))
    write.table(atlas, outfile, sep = ",",
                quote = FALSE, col.names = NA)
  }

  return(as.matrix(atlas))
}
