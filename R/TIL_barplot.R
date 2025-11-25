#' Barplot of abundance estimates
#'
#' Draw barplot of the "betas" from a decon fit
#'
#' @param mat Matrix of cell proportions or abundances, in the same dimensions
#' output by spatialdecon
#'  (cells in rows, observations in columns). User is free to re-order
#'  columns/observations in
#'  whatever order is best for display.
#' @param draw_legend Logical. If TRUE, the function draws a legend in a new
#' plot frame.
#' @param main Title for barplot
#' @param col Vector of colors for cell types. Defaults to pre-set colors for
#' the safeTME cell types.
#' @param ... Arguments passed to barplot()
#' @return Draws a barplot.
#' @examples
#' data(mini_geomx_dataset)
#' data(safeTME)
#' # estimate background:
#' mini_geomx_dataset$bg <- derive_GeoMx_background(
#'   norm = mini_geomx_dataset$normalized,
#'   probepool = rep(1, nrow(mini_geomx_dataset$normalized)),
#'   negnames = "NegProbe"
#' )
#' # run basic decon:
#' res0 <- bulkdecon(
#'   norm = mini_geomx_dataset$normalized,
#'   bg = mini_geomx_dataset$bg,
#'   X = safeTME
#' )
#' # run barplot:
#' TIL_barplot(mat = res0$beta)
#' # run barplot and draw a color legend
#' TIL_barplot(mat = res0$beta, draw_legend = TRUE)
#' @import graphics
#' @importFrom graphics frame legend barplot
#' @importFrom grDevices colors
#' @importFrom utils data
#' @export
TIL_barplot <- function(mat, draw_legend = FALSE, main = "", col = NULL, ...) {


    # infer colors:
    if (length(col) == 0) {

        # use safeTME colors if the right cells are present:
        #utils::data("cellcols", envir = environment())
        # if (all(is.element(rownames(mat), names(SpatialDecon::cellcols)))) {
        #     col <- SpatialDecon::cellcols[rownames(mat)]
        # }
        # else {
            manycols <- c(
                "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
                "#B3DE69", "#FCCDE5", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
                "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#1B9E77", "#D95F02",
                "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
                sample(grDevices::colors(), 99)
            )
            col <- manycols[seq_len(nrow(mat))]
            names(col) <- rownames(mat)
        # }
    }

    usecells <- rownames(mat)

    # draw barplot:
    graphics::barplot(mat[usecells, ],
                      cex.lab = 1.5,
                      col = col, border = NA,
                      las = 2, main = main, ...
    )

    # draw a legend:
    if (draw_legend) {
        graphics::frame()
        graphics::legend("center",
                         fill = rev(col),
                         legend = rev(names(col))
        )
    }
}
