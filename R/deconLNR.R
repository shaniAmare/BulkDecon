#' Log-normal regression deconvolution for bulk RNA-seq
#'
#' Performs deconvolution of bulk RNA-seq expression using a log-normal
#' regression model. This version removes all NanoString-specific assumptions
#' and is designed for bulk RNA-seq abundance matrices only.
#'
#' @param Y Numeric matrix or vector of bulk gene expression values (genes × samples).
#' @param X Numeric matrix of reference gene expression profiles
#'           (genes × cell types), typically derived from single-cell data.
#' @param bg Optional background matrix of same dimension as Y, or a single
#'           scalar (default = 0).
#' @param weights Optional matrix of per-gene weights (inverse variance).
#'                If NULL, equal weights are used.
#' @param epsilon Small positive constant used to stabilize log-normal fitting.
#'                If NULL, epsilon is automatically selected per sample.
#' @param maxit Maximum iterations for the optimizer (default = 1000).
#'
#' @return A list containing:
#'   \item{beta}{Matrix of estimated cell-type abundances (cell types × samples)}
#'   \item{sigmas}{Array of covariance matrices for each sample}
#'
#' @export
deconLNR <- function(Y, X, bg = 0, weights = NULL, epsilon = NULL, maxit = 1000) {

  # -------------------------------
  # Input formatting
  # -------------------------------
  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
    colnames(Y) <- "Sample1"
  }
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.matrix(X)) X <- as.matrix(X)

  if (nrow(Y) != nrow(X)) {
    stop("Y and X must have the same number of genes (rows).")
  }

  # background handling
  if (length(bg) == 1) {
    bg <- matrix(bg, nrow(Y), ncol(Y))
  }

  # weights
  if (is.null(weights)) {
    weights <- matrix(1, nrow(Y), ncol(Y))
  }

  # choose epsilon
  if (is.null(epsilon)) {
    epsilon <- apply(Y, 2, function(y) {
      min(y[y > 0], na.rm = TRUE)
    })
  }

  # -------------------------------
  # Function applied per sample
  # -------------------------------
  decon_single <- function(y, b, wts, eps) {

    # remove NA rows
    keep <- !is.na(y)
    y <- y[keep]
    b <- b[keep]
    wts <- wts[keep]
    Xs <- X[keep, , drop = FALSE]

    # initial values
    init <- rep(mean(y) / (mean(Xs) * ncol(Xs)), ncol(Xs))
    names(init) <- colnames(Xs)

    # fit log-normal model
    fit <- logNormReg::lognlm(
      pmax(y, eps) ~ b + Xs - 1,
      lik = FALSE,
      weights = wts,
      start = c(1, init),
      method = "L-BFGS-B",
      lower = c(1, rep(0, ncol(Xs))),
      upper = c(10, rep(Inf, ncol(Xs))),
      opt = "optim",
      control = list(maxit = maxit)
    )

    list(
      beta = fit$coefficients[-1],
      sigma = solve(fit$hessian)[-1, -1]
    )
  }

  # -------------------------------
  # Apply across samples
  # -------------------------------
  fnlist <- lapply(seq_len(ncol(Y)), function(i) {
    decon_single(
      y = Y[, i],
      b = bg[, i],
      wts = weights[, i],
      eps = epsilon[i]
    )
  })

  # extract beta and sigma
  beta <- sapply(fnlist, function(z) z$beta)
  rownames(beta) <- colnames(X)
  colnames(beta) <- colnames(Y)

  sigmas <- array(
    data = unlist(lapply(fnlist, function(z) z$sigma)),
    dim = c(ncol(X), ncol(X), ncol(Y)),
    dimnames = list(colnames(X), colnames(X), colnames(Y))
  )

  list(
    beta = pmax(beta, 0),
    sigmas = sigmas
  )
}
