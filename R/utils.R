#' @keywords internal
.auto_colors <- function(n, names = NULL) {
  base_cols <- c(
    "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
    "#B3DE69", "#FCCDE5", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
    "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#1B9E77", "#D95F02",
    "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
    "darkred", "violet", "darkblue", "#654522", "#dcd300", "blue",
    "darkgreen"
  )

  # extend if needed
  if (n > length(base_cols)) {
    extra <- grDevices::colors()[grep("^#", grDevices::colors(), invert = TRUE)]
    base_cols <- c(base_cols, sample(extra, n - length(base_cols)))
  }

  cols <- base_cols[seq_len(n)]
  if (!is.null(names)) names(cols) <- names
  cols
}
