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
#'   \item \code{calc_bgindex}: Internal function to compute background-stability
#'         indices (SEG-like) via Normal/Gamma mixture modelling.
#'   \item \code{calc_background}: Uses `calc_bgindex` to identify robust background
#'         genes and estimate sample-wise background levels from a normalized matrix.
#' }
#'
#' @section Functions to set up deconvolution:
#' \itemize{
#'   \item \code{collapseCellTypes}: Merge granular cell types into broader groups.
#'   \item \code{download_profile_matrix}: Download or load profile matrices.
#'   \item \code{safeTME}: Dataset providing immune cell expression profiles for
#'         tumor/immune deconvolution.
#' }
#'
#' @section Deconvolution functions:
#' \itemize{
#'   \item \code{bulkdecon}: Core function performing bulk RNA-seq deconvolution
#'         using supplied cell-type signatures.
#'   \item \code{reverseDecon}: Fits a reverse model (expression ~ cell
#'         proportions) to inspect gene-level dependence on cellular composition
#'         and compute de-mixed residuals.
#' }
#'
#' @section Plotting functions:
#' \itemize{
#'   \item \code{florets}: Visualise cell proportions using cockscomb-style plots
#'         positioned in 2D space.
#'   \item \code{TIL_barplot}: Plot estimated tumor infiltrating lymphocyte (TIL)
#'         abundances from a SafeTME-derived deconvolution.
#' }
#'
#' @docType package
#' @name BulkDecon
NULL
