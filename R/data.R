#' Mapping from granular cell populations to broader cell populations
#'
#' This dataset defines how fine-grained single-cell populations map to
#' broader immune or stromal categories. It is used internally by
#' `convertCellTypes()` to aggregate cell types.
#'
#' @format A named list. Each list element corresponds to a broad cell
#'   population, whose value is a character vector of granular sub-populations.
#'
#' @docType data
#' @usage NULL
#' @keywords datasets
"safeTME.matches"

#' SafeTME Reference Matrix
#'
#' Reference expression matrix of 906 genes across 18 immune and stromal
#' cell types, adapted from the SafeTME reference panel.
#'
#' @format A numeric matrix with 906 rows (genes)
#'   and 18 columns (cell types).
#'
#' @docType data
#' @usage NULL
#' @keywords datasets
"safeTME"

#' Mini Human Colon Single-Cell Dataset
#'
#' A small single-cell RNA-seq dataset containing 250 randomly selected
#' cells from the Kinchen et al. (2018) human colon study. Only genes with
#' CV > 10 across cell types were retained.
#'
#' @format A list with two elements:
#'   \itemize{
#'     \item \code{mtx}: sparse count matrix (genes × cells)
#'     \item \code{annots}: data frame of cell-type annotations
#'   }
#'
#' @docType data
#' @usage NULL
#' @keywords datasets
"mini_singleCell_dataset"

#' Default colors for the cell types in the safeTME matrix
#'
#' A named vector of colors, giving colors for the cell types of the safeTME
#'  matrix.
#'
#' @format A named vector
"cellcols"

#' Small example GeoMx data
#'
#' A miniature GeoMx dataset used by the spatialdecon examples.
#'
#' @format A list with the following elements:
#'  \itemize{
#'  \item normalized: normalized data matrix
#'  \item raw: raw data matrix
#'  \item annot: AOI annotation data frame
#'  }
"mini_geomx_dataset"

#' Genes' biological variability in immune deconvolution from TCGA.
#'
#' Genes' biological SDs, as estimated from immune deconvolution from TCGA.
#' Used to weight genes in spatialdecon.
#'
#' @format A named vector giving SDs of 1179 genes.
"mean.resid.sd"

