#' @title Floret Plot for Compositional Data
#'
#' @description
#' Draws “florets” (radial polygon slices) at given spatial coordinates to
#' visualise compositional vectors (e.g., cell-type proportions) per point.
#' Each row of `b` is drawn as a petal around each `(x, y)` coordinate.
#'
#' @param x Numeric vector of x-coordinates.
#' @param y Numeric vector of y-coordinates.
#' @param b Vector or matrix of values to display as florets.
#'   Rows = categories (e.g., cell types).
#'   Columns = samples matching `x` and `y`.
#' @param col Optional vector or matrix of colours.
#'   If missing, colours are automatically generated.
#' @param legendwindow Logical; if `TRUE`, draws a separate legend-only window of florets.
#' @param rescale.by.sqrt Logical; if `TRUE` (default) the radii are rescaled by the square root.
#' @param border Border colour for floret polygons.
#' @param add Logical; if `TRUE`, adds florets onto an existing plot.
#' @param cex Numeric size scaling factor for florets.
#' @param bty,xlab,ylab,xaxt,yaxt Additional graphical parameters passed to `plot()`.
#' @param ... Additional parameters passed to base graphics.
#'
#' @details
#' This function creates a radial flower-like visualisation, where each sample
#' is represented by a set of coloured polygon slices (“florets”), each slice
#' corresponding to one category (e.g., a deconvolution proportion).
#'
#' Useful for displaying multiple composition vectors at their spatial locations,
#' such as cell-type proportions mapped onto imaging-based spatial coordinates.
#'
#' @return
#' Invisibly returns `NULL`. The function is called for its side-effect of
#' creating a plot.
#'
#' @export
#'
#' @examples
#' ## Example data
#' set.seed(1)
#' x <- runif(5)
#' y <- runif(5)
#' b <- matrix(runif(5 * 3), nrow = 3,
#'             dimnames = list(c("A", "B", "C"), NULL))
#'
#' florets(x, y, b)
florets <- function(x, y, b, col = NULL, legendwindow = FALSE,
                    rescale.by.sqrt = TRUE, border = NA, add = FALSE, cex = 1,
                    bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                    ...) {

  if (rescale.by.sqrt) {
    b <- sqrt(b)
  }

  if (is.vector(b)) {
    b2 <- matrix(b, nrow = length(b))
    rownames(b2) <- names(b)
    b <- b2
    rm(b2)
  }

  # harmonised colour palette
  if (is.null(col) || length(col) == 0) {
    col <- .auto_colors(nrow(b), names = rownames(b))
  }

  if (length(col) == 1)
    col <- matrix(col, nrow = nrow(b), ncol = ncol(b))

  if (is.vector(col))
    col <- matrix(col, nrow = nrow(b), ncol = ncol(b))

  angles <- seq(0, 2 * pi, length.out = nrow(b) + 1)

  maxrange <- max(diff(range(x, na.rm = TRUE)), diff(range(y, na.rm = TRUE)))
  b <- b * maxrange / mean(b, na.rm = TRUE) * 0.007 * cex

  if (!add) {
    graphics::plot(x, y,
                   col = 0, bty = bty, xaxt = xaxt, yaxt = yaxt,
                   xlab = xlab, ylab = ylab, ...)
  }

  if (nrow(b) > 1) {
    for (i in seq_len(length(x))) {
      for (j in seq_len(nrow(b))) {
        tempangles <- seq(angles[j], angles[j + 1], length.out = 20)
        xt <- b[j, i] * cos(tempangles)
        yt <- b[j, i] * sin(tempangles)
        graphics::polygon(
          x[i] + c(0, xt), y[i] + c(0, yt),
          col = col[j, i], border = border, lwd = 0.5
        )
      }
    }
  }

  if (nrow(b) == 1) {
    for (i in seq_len(length(x))) {
      for (j in seq_len(nrow(b))) {
        tempangles <- seq(angles[j], angles[j + 1], length.out = 20)
        xt <- b[j] * cos(tempangles)
        yt <- b[j] * sin(tempangles)
        graphics::polygon(
          x[i] + xt, y[i] + yt,
          col = col[j], border = border, lwd = 0.5
        )
      }
    }
  }

  if (legendwindow) {
    graphics::plot(0, 0,
                   col = 0, xlim = c(-1, 1), ylim = c(-1, 1), xaxt = "n",
                   yaxt = "n", xlab = "", ylab = "", ...
    )
    for (j in seq_len(length(angles))) {
      graphics::lines(
        c(0, 0.75 * cos(angles[j])), c(0, 0.75 * sin(angles[j])),
        col = col[j], lwd = 2
      )
      graphics::text(
        0.85 * cos(angles[j]), 0.85 * sin(angles[j]),
        rownames(b)[j], cex = 1.4
      )
    }
  }
}
