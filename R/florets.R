
florets <- function(x, y, b, col = NULL, legendwindow = FALSE,
                    rescale.by.sqrt = TRUE, border = NA, add = FALSE, cex = 1,
                    bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                    ...) {
    # utils::data("cellcols", envir = environment())

    if (rescale.by.sqrt) {
        b <- sqrt(b)
    }
    # make b a matrix:
    if (is.vector(b)) {
        b2 <- matrix(b, nrow = length(b))
        rownames(b2) <- names(b)
        b <- b2
        rm(b2)
    }
    # choose colors if not given:
    # if ((length(col) == 0) &
    #     all(is.element(rownames(b), names(BulkDecon::cellcols)))) {
    #     col <- BulkDecon::cellcols[rownames(b)]
    # }
    # if ((length(col) == 0) &
    #     !all(is.element(rownames(b), names(BulkDecon::cellcols)))) {
        manycols <- c(
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
          "#B3DE69", "#FCCDE5", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
          "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#1B9E77", "#D95F02",
          "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
          "darkred",  "violet", "darkblue",  "#654522", "#dcd300","blue",
          "darkgreen"
        )
        col <- manycols[seq_len(nrow(b))]
        names(col) <- rownames(b)
    # }
    # convert colors to matrix of the same dimension as b:
    if (length(col) == 1) {
        col <- matrix(col, nrow = nrow(b), ncol = ncol(b))
    }
    if (is.vector(col)) {
        col <- matrix(col, nrow = nrow(b), ncol = ncol(b))
    }

    # get radians:
    angles <- seq(0, 2 * pi, length.out = nrow(b) + 1)

    # scale b based on the range of x and y:
    maxrange <- max(diff(range(x, na.rm = TRUE)), diff(range(y, na.rm = TRUE)))
    b <- b * maxrange / mean(b, na.rm = TRUE) * 0.007 * cex

    # draw plot:
    if (!add) {
        graphics::plot(x, y,
                       col = 0, bty = bty, xaxt = xaxt, yaxt = yaxt,
                       xlab = xlab, ylab = ylab, ...
        )
    }

    # draw florets:
    if (nrow(b) > 1) {
        for (i in seq_len(length(x))) {
            for (j in seq_len(nrow(b))) {
                tempangles <- seq(angles[j], angles[j + 1], length.out = 20)
                xt <- b[j, i] * cos(tempangles)
                yt <- b[j, i] * sin(tempangles)
                graphics::polygon(x[i] + c(0, xt), y[i] + c(0, yt),
                                  col = col[j, i],
                                  border = border, lwd = 0.5
                )
            }
        }
    }

    # if just one point, draw a full circle:
    if (nrow(b) == 1) {
        for (i in seq_len(length(x))) {
            for (j in seq_len(nrow(b))) {
                tempangles <- seq(angles[j], angles[j + 1], length.out = 20)
                xt <- b[j, i] * cos(tempangles)
                yt <- b[j, i] * sin(tempangles)
                graphics::polygon(x[i] + xt, y[i] + yt,
                                  col = col[j],
                                  border = border, lwd = 0.5
                )
            }
        }
    }

    # draw a legend:
    if (legendwindow) {
        graphics::plot(0, 0,
                       col = 0, xlim = c(-1, 1), ylim = c(-1, 1), xaxt = "n",
                       yaxt = "n", xlab = "", ylab = "", ...
        )
        for (j in seq_len(length(angles))) {
            graphics::lines(c(0, 0.75 * cos(angles[j])), c(0, 0.75 * sin(angles[j])),
                            col = col[j], lwd = 2
            )
            graphics::text(0.85 * cos(angles[j]), 0.85 * sin(angles[j]),
                           rownames(b)[j],
                           cex = 1.4
            )
        }
    }
}
