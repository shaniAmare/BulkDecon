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
runErrorModel <- function(counts, method = "quantile") {

  if (!is.numeric(counts)) stop("counts must be numeric.")
  method <- tolower(method)

  # --------------------------------------
  # Quantile-based empirical variance model
  # --------------------------------------
  if (method == "quantile") {

    if (is.vector(counts)) {
      quant <- rank(counts, ties.method = "average") / length(counts)
    } else {
      quant <- apply(counts, 2, function(x) rank(x) / length(x))
    }

    # baseline variance is small
    sds <- quant * 0 + 0.1

    # lower abundance → more noise
    sds[quant < 0.2]  <- 0.2
    sds[quant < 0.15] <- 0.3
    sds[quant < 0.10] <- 0.4
    sds[quant < 0.05] <- 0.5
    sds[quant < 0.01] <- 1.0

    return(sds)
  }

  # --------------------------------------
  # Log-linear abundance–variance model
  # --------------------------------------
  if (method == "loglinear") {

    logc <- log2(pmax(counts, 1))

    # sd shrinks with log abundance
    sds <- 1.5 - 0.2 * logc
    sds[sds < 0.1] <- 0.1

    return(sds)
  }

  stop("Unknown method: choose 'quantile' or 'loglinear'.")
}
