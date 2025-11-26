#' BulkDecon: Bulk RNA-seq Deconvolution Using Custom Single-Cell References
#'
#' The **BulkDecon** package performs bulk RNA-seq deconvolution using user-defined
#' single-cell reference datasets. It provides tools to preprocess expression
#' matrices, estimate background levels, construct cell-type signature matrices,
#' run deconvolution models, inspect residuals, and visualise results.
#'
#' @section Background estimation helpers:
#' Functions for identifying stably low-expressed genes and estimating
#' sample-specific background for bulk RNA-seq or GeoMx-style data.
#'
#' \itemize{
#'   \item \code{calc_bgindex}: Compute background-stability indices via
#'         Normal/Gamma mixture modelling.
#'   \item \code{calc_background}: Estimate robust background levels
#'         across samples.
#' }
#'
#' @section Functions to set up deconvolution:
#' \itemize{
#'   \item \code{collapseCellTypes}: Merge granular cell types.
#'   \item \code{download_profile_matrix}: Download or load profile matrices.
#'   \item \code{safeTME}: Immune profile reference.
#' }
#'
#' @section Deconvolution functions:
#' \itemize{
#'   \item \code{bulkdecon}: Perform bulk RNA-seq deconvolution.
#'   \item \code{reverseDecon}: Fit reverse model and calculate residuals.
#' }
#'
#' @section Plotting functions:
#' \itemize{
#'   \item \code{florets}: Visualise cell proportions.
#'   \item \code{TIL_barplot}: Tumour-infiltrating lymphocyte barplots.
#' }
#'
#' @keywords internal
"_PACKAGE"
