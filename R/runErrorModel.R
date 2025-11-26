#' Estimate gene-specific standard deviations for bulk RNA-seq
#'
#' Produces per-gene SD estimates used for weighting in log-normal regression
#' deconvolution. Uses a quantile-based empirical model by default.
#'
#' @param counts Numeric vector or matrix of raw or normalized bulk RNA-seq counts.
#' @param method Method for error modelling. One of:
#'   - "quantile" (default): empirical stability by rank
#'   - "loglinear": sd increases with decreasing abundance via log model
#'
#' @return A numeric vector or matrix of estimated standard deviations.
#'
#' @export
runErrorModel <- function(counts, platform = "general") {
  if (platform == "ncounter") {
    sds <- counts * 0 + 0.1
    sds <- replace(sds, counts < 200, 0.2)
    sds <- replace(sds, counts < 100, 0.3)
    sds <- replace(sds, counts < 75, 0.4)
    sds <- replace(sds, counts < 50, 0.5)
    sds <- replace(sds, counts < 40, 0.7)
    sds <- replace(sds, counts < 30, 1)
    sds <- replace(sds, counts < 20, 3)
  }


  if (platform == "rsem") {
    sds <- counts * 0 + 0.5930982
    sds <- replace(sds, log2(counts) < 9.5, 0.6458475)
    sds <- replace(sds, log2(counts) < 8.5, 0.7847597)
    sds <- replace(sds, log2(counts) < 7.5, 1.0576471)
    sds <- replace(sds, log2(counts) < 6.5, 1.2990917)
    sds <- replace(sds, log2(counts) < 5.5, 1.5061735)
    sds <- replace(sds, log2(counts) < 4.5, 1.6930872)
    sds <- replace(sds, log2(counts) < 3.5, 1.7894239)
  }

  if (platform == "bulk") {
    predictsd.dsp <- function(rawcounts) {
      m <- log2(pmax(rawcounts, 1e-3))
      meanvec <- c(-1e-6, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, Inf)
      sdvec <- c(
        1.5, 1.383, 1.191, 0.800, 0.48, 0.301, 0.301,
        0.301, 0.263, 0.235, 0.235
      )

      s <- replace(m, TRUE, sdvec[1])
      for (i in seq_len(length(meanvec) - 1)) {
        s <- replace(s, m >= meanvec[i], sdvec[i + 1])
      }
      return(s)
    }

    if (is.vector(counts)) {
      sds <- vapply(
        X = counts,
        FUN = predictsd.dsp,
        FUN.VALUE = numeric(length(counts))
      )
    }
    if (is.matrix(counts)) {
      sds <- predictsd.dsp(counts)
    }
  }

  if (platform == "st") {
    # assume poisson error:
    sds <- counts * 0 + 03
    sds <- replace(sds, counts < 500, 0.045)
    sds <- replace(sds, counts < 200, 0.07)
    sds <- replace(sds, counts < 100, 0.1)
    sds <- replace(sds, counts < 50, 0.14)
    sds <- replace(sds, counts < 30, 0.18)
    sds <- replace(sds, counts < 20, 0.23)
    sds <- replace(sds, counts < 15, 0.27)
    sds <- replace(sds, counts < 10, 0.35)
    sds <- replace(sds, counts < 5, 0.61)
    sds <- replace(sds, counts < 2, 1.15)
    sds <- replace(sds, counts < 1, 1.33)
  }

  if (platform == "quantile") {
    if (is.vector(counts)) {
      quantile <- rank(counts) / length(counts)
    }
    if (is.matrix(counts)) {
      quantile <- matrix(rank(counts), nrow(counts)) / length(counts)
    }

    sds <- quantile * 0 + 0.1
    sds <- replace(sds, quantile < 0.2, 0.2)
    sds <- replace(sds, quantile < 0.15, 0.3)
    sds <- replace(sds, quantile < 0.1, 0.4)
    sds <- replace(sds, quantile < 0.05, 0.5)
    sds <- replace(sds, quantile < 0.01, 1)
  }
  return(sds)
}
