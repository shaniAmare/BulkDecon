#' Plot TIL Composition as Barplots
#'
#' @description
#' Creates a stacked barplot of cell-type proportions across samples.
#' Each row of `mat` represents a cell type, and each column represents a sample.
#'
#' @details
#' If no color vector is provided, the function automatically generates a large
#' set of distinct colors for all cell types.
#' A legend can be optionally drawn beside the barplot.
#'
#' @param mat A numeric matrix with **cell types in rows** and **samples in columns**.
#'   Row names must be present and will be used as labels.
#' @param draw_legend Logical; if `TRUE`, draws a legend in a separate frame.
#' @param main Character; title of the barplot.
#' @param col Optional named vector of colors, one per cell type (rownames of `mat`).
#'   If omitted, a reproducible palette is auto-generated.
#' @param ... Additional arguments passed to `graphics::barplot()`.
#'
#' @return
#' Invisibly returns `NULL`.
#' The function is used for generating plots and has no structured return object.
#'
#' @examples
#' \dontrun{
#' mat <- matrix(runif(30), nrow = 5)
#' rownames(mat) <- c("T_cells", "B_cells", "NK", "Macrophages", "DCs")
#' colnames(mat) <- paste0("Sample", 1:6)
#'
#' TIL_barplot(mat, draw_legend = TRUE, main = "TIL Composition")
#' }
#'
#' @export
TIL_barplot <- function(mat, draw_legend = FALSE, main = "", col = NULL, ...) {

  if (!is.matrix(mat))
    stop("`mat` must be a matrix with cell types in rows and samples in columns.")

  if (is.null(rownames(mat)))
    stop("`mat` must have rownames corresponding to cell types.")

  # harmonised colours
  if (is.null(col) || length(col) == 0) {
    col <- .auto_colors(nrow(mat), names = rownames(mat))
  }

  usecells <- rownames(mat)

  graphics::barplot(
    mat[usecells, , drop = FALSE],
    cex.lab = 1.5,
    col = col[usecells],
    border = NA,
    las = 2,
    main = main,
    ...
  )

  if (draw_legend) {
    graphics::frame()
    graphics::legend(
      "center",
      fill = rev(col[usecells]),
      legend = rev(usecells)
    )
  }
}
