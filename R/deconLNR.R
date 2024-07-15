
deconLNR <- function(Y, X, bg = 0, weights = NULL, epsilon = NULL,
                     maxit = 1000) {
    if (length(weights) == 0) {
        weights <- replace(Y, TRUE, 1)
    }
    if (is.vector(Y)) {
        Y <- matrix(Y, nrow = length(Y))
    }
    if (length(bg) == 1) {
        bg <- matrix(bg, nrow(Y), ncol(Y))
    }
    # choose "epsilon": a very small non-zero number to make fits well-behaved
    if (length(epsilon) == 0) {
        epsilon <- min(replace(Y, (Y == 0) & !is.na(Y), NA), na.rm = TRUE)
    }
    
    # matrix-like data for apply:
    mat <- rbind(Y, bg, weights)
    # fn to apply:
    fn <- function(zz) {
        # break into y, b, w:
        y <- zz[seq_len((length(zz)) / 3)]
        b <- zz[seq_len((length(zz)) / 3) + (length(zz) / 3)]
        wts <- zz[seq_len((length(zz)) / 3) + (length(zz) / 3) * 2]
        
        # remove NA data:
        use <- !is.na(y)
        y <- y[use]
        b <- b[use]
        Xtemp <- X[use, , drop = FALSE]
        wts <- wts[use]
        
        init <- rep(mean(y) / (mean(X) * ncol(X)), ncol(X))
        names(init) <- colnames(X)
        
        # run lognlm:
        fit <- logNormReg::lognlm(pmax(y, epsilon) ~ b + Xtemp - 1,
                       lik = FALSE,
                       weights = wts,
                       start = c(1, init),
                       method = "L-BFGS-B",
                       lower = c(1, rep(0, ncol(Xtemp))),
                       upper = c(1, rep(Inf, ncol(Xtemp))),
                       opt = "optim",
                       control = list(maxit = maxit)
        )
        fnout <- list(
            beta = fit$coefficients[-1],
            sigma = solve(fit$hessian)[-1, -1]
        )
        return(fnout)
    }
    # apply across all observations:
    fnlist <- apply(mat, 2, fn)
    # extract beta and sigmas:
    getbeta <- function(zz) {
        return(zz$beta)
    }
    getsigma <- function(zz) {
        return(zz$sigma)
    }
    
    beta <- vapply(X = fnlist, FUN = getbeta, FUN.VALUE = numeric(ncol(X)))
    rownames(beta) <- colnames(X)
    sigmas <- array(vapply(
        X = fnlist,
        FUN = getsigma,
        FUN.VALUE = numeric(ncol(X)^2)
    ),
    dim = c(ncol(X), ncol(X), ncol(Y)),
    dimnames = list(colnames(X), colnames(X), colnames(Y))
    )
    
    out <- list(beta = pmax(beta, 0), sigmas = sigmas)
    return(out)
}
