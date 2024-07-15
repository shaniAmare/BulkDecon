#' @section functions:
#' Functions to help set up deconvolution:
#' \itemize{
#'  \item derive_GeoMx_background Estimates the background levels from GeoMx
#'  experiments
#'  \item collapseCellTypes reformats deconvolution results to merge
#'  closely-related cell types
#'  \item download_profile_matrix Downloads a cell profile matrix.
#'  \item safeTME: a data object, a matrix of immune cell profiles for use in
#'   tumor-immune deconvolution.
#' }
#' Deconvolution functions:
#' \itemize{
#'  \item spatialdecon runs the core deconvolution function
#'  \item reverseDecon runs a transposed/reverse deconvolution problem, fitting
#'  the data as a function of cell abundance estimates.
#'   Used to measure genes' dependency on cell mixing and to calculate gene
#'    residuals from cell mixing.
#' }
#' Plotting functions:
#' \itemize{
#'  \item florets Plot cell abundance on a specified x-y space, with each point
#'   a cockscomb plot showing the cell abundances of that region/sample.
#'  \item TIL_barplot Plot abundances of tumor infiltrating lymphocytes (TILs)
#'   estimated from the safeTME cell profile matrix
#' }
#' @docType package
#' @name BulkDecon-package
NULL
