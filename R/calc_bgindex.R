#' Calculate background-stability index for genes
#'
#' Internal function used to identify genes that are stably low (background)
#' across samples. The function fits a two-component Normal/Gamma mixture model
#' (EM) to non-zero expression values per gene and combines three stability
#' features into a single `bgIdx` score.
#'
#' @param exprs_mat Numeric matrix (genes × samples). Values must be >= 0.
#' @param BPPARAM BiocParallel parameter object; default uses SerialParam.
#' @param return_all Logical; if TRUE, returns extra columns including zeros for
#'   all input genes (not only filtered genes).
#'
#' @return If \code{return_all = FALSE} (default) a data.frame with rows for
#'   genes that passed the zero-filter and columns: bgIdx, rho, sigma, mu,
#'   mu.scaled, zero. If \code{return_all = TRUE} returns a data.frame with one
#'   row per input gene including zero fraction.
#'
#' @details
#' Steps:
#' 1. Compute fraction of zeros per gene and remove genes with >80% zeros.
#' 2. For remaining genes, fit a 2-component Normal/Gamma mixture (zeroes
#'    removed) using \code{gammaNormMix} in parallel.
#' 3. Combine scaled mixture parameters (rho, sigma, mu) and zero fraction
#'    into a composite background index \code{bgIdx} (higher = more background-like).
#'
#' @keywords internal
#' @importFrom BiocParallel SerialParam bplapply
#' @export
calc_bgindex <- function(exprs_mat,
                         BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
                         return_all = FALSE) {

  ## Input checks
  if (is.null(exprs_mat)) stop("`exprs_mat` is NULL.")
  if (!is.matrix(exprs_mat)) stop("`exprs_mat` must be a matrix (genes × samples).")
  if (any(exprs_mat < 0, na.rm = TRUE)) stop("`exprs_mat` contains negative values.")

  # compute zero fraction BEFORE any log transform
  zero_frac_all <- rowMeans(exprs_mat == 0, na.rm = TRUE)

  # filter genes that are > 80% zero
  del <- which(zero_frac_all > 0.8)
  if (length(del) > 0) {
    exprs_mat_filt <- exprs_mat[-del, , drop = FALSE]
    if (nrow(exprs_mat_filt) < 100) {
      stop("Not enough genes pass QC after zero-filtering (< 100).")
    }
  } else {
    exprs_mat_filt <- exprs_mat
  }

  message("Calculating background index...")

  # take log of positive values for modelling (preserve zeros removed above)
  log_exprs <- log(exprs_mat_filt)
  # replace -Inf by 0? We removed >80% zeros genes, but zeros remain; gammaNormMix handles removeZeroes
  # call parallel parameter to compute mixture parameters per gene
  paraMat <- make_para_gn_parallel(as.matrix(log_exprs), BPPARAM = BPPARAM)
  r <- paraMat$rho
  s <- paraMat$sigma
  m <- paraMat$mu

  # scale mu to 0..1
  m.scaled <- (m - min(m, na.rm = TRUE)) / (max(m, na.rm = TRUE) - min(m, na.rm = TRUE))
  names(r) <- names(s) <- names(m) <- names(m.scaled) <- rownames(exprs_mat_filt)

  genes <- rownames(exprs_mat_filt)
  # combine zero fraction (from original data) with scaled mu
  z <- zero_frac_all[genes] * m.scaled[genes]

  # features to combine (rank-based)
  x1 <- rank(r, ties.method = "average") / (length(r) + 1)
  x2 <- 1 - rank(s, ties.method = "average") / (length(s) + 1)
  x3 <- 1 - rank(z, ties.method = "average") / (length(z) + 1)

  bgIdx <- rowMeans(cbind(x1, x2, x3), na.rm = TRUE)

  resMat <- data.frame(bgIdx = bgIdx, rho = r, sigma = s,
                       mu = m, mu.scaled = m.scaled, zero = z,
                       stringsAsFactors = FALSE)
  rownames(resMat) <- genes

  if (return_all) {
    resMat2 <- cbind(gene = as.character(rownames(resMat)), resMat)
    all_genes <- data.frame(gene = as.character(rownames(exprs_mat)),
                            zero.all = zero_frac_all,
                            stringsAsFactors = FALSE)
    result <- merge(x = all_genes, y = resMat2, by = "gene", all.x = TRUE)
    rownames(result) <- rownames(exprs_mat)
  } else {
    result <- resMat
  }

  return(result)
}


#' Internal helper: compute mixture parameters for a single gene
#'
#' @keywords internal
compute_data_para <- function(exprs_mat_gene) {
  # exprs_mat_gene: numeric vector (log-expr for a single gene)
  # run gammaNormMix on non-NA values
  if (all(is.na(exprs_mat_gene))) {
    return(c(mu = NA, sigma = NA, rho = NA))
  }
  para <- tryCatch(
    gammaNormMix(exprs_mat_gene, verbose = FALSE, plot = FALSE, removeZeroes = TRUE, onlyAddCurves = TRUE),
    error = function(e) NA
  )

  if (is.list(para) && !is.null(para$mu) && !is.na(para$mu)) {
    data_para <- c(mu = para$mu, sigma = para$sd, rho = para$rho)
  } else if (is.list(para) && !is.null(para$mu) && !is.na(para$mu)) {
    data_para <- c(mu = para$mu, sigma = para$sd, rho = para$rho)
  } else {
    # fallback NA vector
    data_para <- c(mu = NA, sigma = NA, rho = NA)
  }
  return(data_para)
}


#' Internal helper: parallel wrapper for compute_data_para
#'
#' @keywords internal
#' @importFrom BiocParallel bplapply
make_para_gn_parallel <- function(exprs_mat, BPPARAM) {
  # exprs_mat: a numeric matrix with genes in rows
  res_list <- BiocParallel::bplapply(
    X = seq_len(nrow(exprs_mat)),
    FUN = function(i) compute_data_para(exprs_mat[i, ]),
    BPPARAM = BPPARAM
  )

  res <- do.call(rbind, res_list)
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  rownames(res) <- rownames(exprs_mat)
  colnames(res) <- c("mu", "sigma", "rho")
  return(res)
}


#' BIC helper
#' @keywords internal
bic <- function(loglik, n, p) {
  -2 * loglik + p * log(n)
}

#' AIC helper
#' @keywords internal
aic <- function(loglik, p) {
  -2 * loglik + 2 * p
}

#' ICL-BIC helper
#' @keywords internal
icl_bic <- function(loglik, postprob, n, p) {
  postprob <- postprob[postprob > 0]
  EN <- -sum(postprob * log(postprob))
  -2 * loglik + 2 * EN + p * log(n)
}


#' Fit a 2-component Normal/Gamma mixture to a vector of values (EM)
#'
#' This function fits a mixture model to non-zero values and returns a list
#' with probExpressed (per-observation posterior) and model parameters.
#' It attempts to reproduce the original behaviour but is more defensive and
#' will fall back gracefully where packages are unavailable.
#'
#' @param data numeric vector
#' @param thresh convergence threshold
#' @param maxiter maximum iterations
#' @param removeZeroes whether to remove zero values before EM
#' @param plot logical; whether to plot diagnostics (default FALSE in compute wrapper)
#' @param verbose logical
#' @return list with elements probExpressed, propExpressed, numExpressed, mu, sd, alpha, beta, rho, niter, loglik, BIC, AIC, ICL_BIC, AreaDifference
#' @keywords internal
gammaNormMix <- function(data,
                         thresh = 1e-07,
                         maxiter = 10000,
                         removeZeroes = TRUE,
                         plot = TRUE,
                         hist = TRUE,
                         hist_col = "light cyan",
                         verbose = FALSE,
                         forceExponential = FALSE,
                         calculateAreaDifference = FALSE,
                         minDataPoints = 5,
                         onlyAddCurves = FALSE,
                         addContextData = FALSE,
                         contextData = NULL) {

  if (addContextData) {
    nOriginal <- length(data)
    data <- c(data, contextData)
  }

  if (removeZeroes) {
    nonZeroInd <- which(data > 0)
    x <- data[nonZeroInd]
  } else {
    x <- data
  }

  if (length(x) < minDataPoints) {
    if (verbose) cat("Not enough data points to fit mixture model!\n")
    return(NA)
  }

  # initialize EM
  n <- length(x)
  z_iter <- stats::rbinom(n, 1, 0.5)
  if (sum(z_iter) == 0) z_iter[1] <- 1

  mu_iter <- mean(x[z_iter == 1], na.rm = TRUE)
  if (is.na(mu_iter)) mu_iter <- mean(x, na.rm = TRUE)
  sig2_iter <- var(x[z_iter == 1], na.rm = TRUE)
  if (is.na(sig2_iter) || sig2_iter <= 0) sig2_iter <- 1e-6
  # initialize gamma params by method of moments on (1 - z) portion or entire data
  alpha_iter <- 2
  beta_iter <- 1
  rho_iter <- mean(z_iter, na.rm = TRUE)

  # safe iterations
  niter <- 0
  converged <- FALSE

  while (niter < maxiter) {
    niter <- niter + 1

    # M-step
    mu_new <- sum(z_iter * x, na.rm = TRUE) / sum(z_iter, na.rm = TRUE)
    sig2_new <- sum(z_iter * (x - mu_new)^2, na.rm = TRUE) / sum(z_iter, na.rm = TRUE)
    if (is.na(sig2_new) || sig2_new <= 0) sig2_new <- 1e-11

    # gamma params: update beta by method of moments approx on (1 - z)
    denom <- sum((1 - z_iter) * x, na.rm = TRUE)
    if (denom <= 0 || is.na(denom)) denom <- sum(x, na.rm = TRUE) + 1e-8
    beta_new <- max(1e-8, alpha_iter * sum(1 - z_iter, na.rm = TRUE) / denom)

    if (!forceExponential) {
      # try to estimate alpha using distr::igamma if available; fallback to 2
      alpha_new <- tryCatch({
        if (requireNamespace("distr", quietly = TRUE)) {
          distr::igamma(sum((log(beta_new) + log(x)) * (1 - z_iter), na.rm = TRUE) / sum(1 - z_iter, na.rm = TRUE))
        } else {
          # fallback: keep alpha at previous value (stable but crude)
          alpha_iter
        }
      }, error = function(e) alpha_iter)
      # cap
      if (is.na(alpha_new) || alpha_new > 150) alpha_new <- min(150, max(0.1, alpha_iter))
    } else {
      alpha_new <- 1
    }

    rho_new <- sum(z_iter, na.rm = TRUE) / n

    # E-step: compute posterior probability z_i = P(component = normal | x_i)
    # logit form:
    eta_iter <- -0.5 * log(2 * pi * sig2_new) - ((x - mu_new)^2) / (2 * sig2_new) -
      alpha_new * log(beta_new) + log(gamma(alpha_new)) - (alpha_new - 1) * log(x + 1e-12) +
      beta_new * x + log(rho_new / (1 - rho_new))

    # avoid overflow: cap eta
    eta_iter[is.infinite(eta_iter)] <- sign(eta_iter[is.infinite(eta_iter)]) * 1e6
    z_iter <- 1 / (1 + exp(-eta_iter))

    # check convergence (parameter-wise)
    if (all(c(abs(mu_new - mu_iter) < thresh,
              abs(sig2_new - sig2_iter) < thresh,
              abs(alpha_new - alpha_iter) < thresh,
              abs(beta_new - beta_iter) < thresh,
              abs(rho_new - rho_iter) < thresh), na.rm = TRUE)) {
      converged <- TRUE
      mu_iter <- mu_new; sig2_iter <- sig2_new; alpha_iter <- alpha_new; beta_iter <- beta_new; rho_iter <- rho_new
      break
    }

    # update params for next iter
    mu_iter <- mu_new; sig2_iter <- sig2_new; alpha_iter <- alpha_new; beta_iter <- beta_new; rho_iter <- rho_new
  }

  # compute log-likelihood
  ll <- sum(log(rho_iter * stats::dnorm(x, mu_iter, sqrt(sig2_iter)) +
                  (1 - rho_iter) * stats::dgamma(x, shape = alpha_iter, rate = beta_iter)),
            na.rm = TRUE)

  # generate smooth curves for plotting if requested
  xg <- seq(0, max(x, na.rm = TRUE) + 1, length.out = 300)
  c1g <- rho_iter * stats::dnorm(xg, mu_iter, sqrt(sig2_iter))
  c2g <- (1 - rho_iter) * stats::dgamma(xg, shape = alpha_iter, rate = beta_iter)
  fg <- c1g + c2g

  if (plot) {
    if (hist) {
      hist(x, probability = TRUE, col = hist_col, breaks = 50,
           main = NA, xlab = NA, ylab = "Density (zeroes removed)",
           ylim = c(0, max(c(c1g, c2g, fg, density(x)$y), na.rm = TRUE)), xlim = c(0, max(x, na.rm = TRUE) + 1))
    }
    if (!onlyAddCurves) {
      graphics::lines(stats::density(x, from = 0), lty = 2, lwd = 2, col = scales::alpha("darkgrey", 0.6))
    }
    graphics::lines(xg, c1g, col = scales::alpha("red", 0.6), lwd = 2)
    graphics::lines(xg, c2g, col = scales::alpha("blue", 0.6), lwd = 2)
    graphics::lines(xg, fg, col = scales::alpha("black", 0.6), lwd = 2)

    if (onlyAddCurves) return(list(xg = xg, c1g = c1g, c2g = c2g, fg = fg))
  }

  if (calculateAreaDifference) {
    f1_fun <- tryCatch({
      dens_fun <- stats::density(x, from = 0)
      dens_interp <- stats::approxfun(dens_fun$x, dens_fun$y)(xg)
      stats::approxfun(xg, dens_interp - fg)
    }, error = function(e) NULL)

    if (!is.null(f1_fun)) {
      AreaDifference <- tryCatch({
        f2 <- function(xx) abs(f1_fun(xx))
        stats::integrate(f2, min(x[x != 0], na.rm = TRUE), max(x, na.rm = TRUE))$value
      }, error = function(e) NULL)
    } else {
      AreaDifference <- NULL
    }
  } else {
    AreaDifference <- NULL
  }

  # reconstruct z for full data length if removed zeros
  if (removeZeroes) {
    z_full <- rep(0, length(data))
    if (exists("nonZeroInd")) z_full[nonZeroInd] <- z_iter
    z <- z_full
  } else {
    z <- z_iter
  }

  # cap probabilities for extremely large values
  maxdata <- data[which.max(z)]
  if (length(maxdata) > 0) {
    z[which(data > maxdata)] <- max(z, na.rm = TRUE)
  }

  # if addContextData, truncate back
  if (addContextData && exists("nOriginal")) {
    z <- z[seq_len(nOriginal)]
  }

  model_bic <- bic(ll, n, 5)
  model_aic <- aic(ll, 5)
  model_icl_bic <- icl_bic(ll, z, n, 5)

  return(list(
    probExpressed = z,
    propExpressed = n * rho_iter / length(data),
    numExpressed = length(which(z > 0.5)),
    mu = mu_iter,
    sd = sqrt(sig2_iter),
    alpha = alpha_iter,
    beta = beta_iter,
    rho = rho_iter,
    niter = niter,
    loglik = ll,
    BIC = model_bic,
    AIC = model_aic,
    ICL_BIC = model_icl_bic,
    AreaDifference = AreaDifference
  ))
}
