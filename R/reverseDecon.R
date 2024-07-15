
reverseDecon <- function(norm, beta, epsilon = NULL) {
  
  # remove cell types with no SD:
  beta = beta[apply(beta, 1, stats::sd) > 0, ]
  # remove cell types that get removed by lm() (presumably removed due to linear dependence)
  lm1 = stats::lm(norm[1,] ~ t(beta))
  beta = beta[!is.na(lm1$coef[setdiff(names(lm1$coef), "(Intercept)")]), , drop = FALSE]
  
  # run reverse decon for all genes:
  rd <- function(y) {
    fit <- suppressWarnings(
      logNormReg::lognlm(y ~ t(beta),
                         lik = FALSE,
                         method = "L-BFGS-B",
                         lower = rep(0, ncol(beta) + 1),
                         upper = rep(Inf, ncol(beta) + 1),
                         opt = "optim",
                         control = list(maxit = 1000)
      )
    )
    return(fit$coefficients)
  }
  coefs <- t(apply(norm, 1, rd))
  colnames(coefs)[-1] <- rownames(beta)
  
  # get yhat
  yhat <- norm * NA
  for (ind in seq_len(ncol(yhat))) {
    yhat[, ind] <- coefs %*% c(1, beta[, ind])
  }
  
  # auto-select a reasonable epsilon if not provided
  if (length(epsilon) == 0) {
    epsilon <- stats::quantile(norm[norm > 0], 0.01)
  }
  
  # get resids:
  resids <- log2(pmax(norm, epsilon)) - log2(pmax(yhat, epsilon))
  
  # get summary stats:
  cors <- suppressWarnings(diag(stats::cor(t(norm), t(yhat))))
  resid.sd <- apply(resids, 1, stats::sd)
  
  out <- list(coefs = coefs, yhat = yhat, resids = resids, cors = cors, resid.sd = resid.sd)
  return(out)
}
