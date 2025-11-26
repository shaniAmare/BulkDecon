#' Reverse Deconvolution of Bulk Gene Expression
#'
#' @description
#' Performs reverse deconvolution of a normalized bulk gene expression matrix
#' given a cell-type reference profile matrix (`beta`).
#' The method fits log-normal linear models (`logNormReg::lognlm`) for each gene,
#' estimates per-cell-type contributions, reconstructs predicted bulk expression,
#' and computes residual diagnostics.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Removes cell types with zero variance across samples.
#'   \item Removes cell types that would be dropped due to linear dependence
#'         (detected via an initial `lm()` call).
#'   \item Fits log-normal regression using \pkg{logNormReg} for each gene.
#'   \item Computes predicted bulk expression \eqn{\hat{Y}}.
#'   \item Selects a reasonable \code{epsilon} for log-transform stability if not supplied.
#'   \item Returns model coefficients, fitted values, residuals, correlations, and residual SD.
#' }
#'
#' @param norm A numeric matrix of normalized gene expression
#'   (genes × samples). Must be non-negative.
#' @param beta A numeric matrix of reference cell-type profiles
#'   (cell types × samples). Row names must match the CTs to estimate.
#' @param epsilon Optional numeric. Lower bound applied before log-transform.
#'   If `NULL`, automatically set to the 1st percentile of positive values in `norm`.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{coefs}{Matrix of model coefficients (intercept + CTs) per gene.}
#'   \item{yhat}{Matrix of predicted gene expression.}
#'   \item{resids}{Matrix of log-residuals.}
#'   \item{cors}{Vector of per-gene correlations between observed and predicted expression.}
#'   \item{resid.sd}{Vector of per-gene residual standard deviations.}
#' }
#'
#' @examples
#' \dontrun{
#' # norm:   gene × sample matrix
#' # beta:   CT × sample matrix
#'
#' result <- reverseDecon(norm = norm, beta = beta)
#'
#' head(result$coefs)
#' plot(result$cors)
#' }
#'
#' @importFrom stats lm sd quantile cor
#' @importFrom logNormReg lognlm
#'
#' @export
reverseDecon <- function(norm, beta, epsilon = NULL) {

  # 1. remove CTs with no SD
  beta <- beta[apply(beta, 1, stats::sd) > 0, , drop = FALSE]

  # 2. remove CTs dropped by lm()
  lm1 <- stats::lm(norm[1, ] ~ t(beta))

  kept <- names(lm1$coef)
  kept <- kept[kept != "(Intercept)"]
  kept <- kept[!is.na(lm1$coef[kept])]

  beta <- beta[kept, , drop = FALSE]

  # 3. Reverse decon per gene
  rd <- function(y) {
    fit <- suppressWarnings(
      logNormReg::lognlm(
        y ~ t(beta),
        lik = FALSE,
        method = "L-BFGS-B",
        lower = rep(0, ncol(beta) + 1),
        upper = rep(Inf, ncol(beta) + 1),
        opt = "optim",
        control = list(maxit = 1000)
      )
    )
    fit$coefficients
  }

  coefs <- t(apply(norm, 1, rd))
  colnames(coefs)[-1] <- rownames(beta)

  # 4. predicted values yhat
  yhat <- matrix(NA, nrow(norm), ncol(norm))
  rownames(yhat) <- rownames(norm)
  colnames(yhat) <- colnames(norm)

  for (i in seq_len(ncol(yhat))) {
    yhat[, i] <- coefs %*% c(1, beta[, i])
  }

  # 5. epsilon selection
  if (is.null(epsilon) || length(epsilon) == 0) {
    positive_vals <- norm[norm > 0]
    epsilon <- if (length(positive_vals) == 0) 1e-6
    else stats::quantile(positive_vals, 0.01)
  }

  # 6. residuals + stats
  resids <- log2(pmax(norm, epsilon)) - log2(pmax(yhat, epsilon))

  cors <- suppressWarnings(diag(stats::cor(t(norm), t(yhat))))
  resid.sd <- apply(resids, 1, stats::sd)

  # 7. return
  list(
    coefs = coefs,
    yhat = yhat,
    resids = resids,
    cors = cors,
    resid.sd = resid.sd
  )
}
