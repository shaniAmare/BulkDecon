
algorithm2 <- function(Y, X, bg = 0, weights = NULL,
                       resid_thresh = 3, lower_thresh = 0.5,
                       align_genes = TRUE, maxit = 1000) {
    
    # align genes:
    if (align_genes) {
        sharedgenes <- intersect(rownames(X), rownames(Y))
        Y <- Y[sharedgenes, ]
        X <- X[sharedgenes, ]
        if (is.matrix(bg)) {
            bg <- bg[sharedgenes, ]
        }
        if (is.matrix(weights)) {
            weights <- weights[sharedgenes, ]
        }
    }
    
    # format the data nicely:
    tidied <- tidy_X_and_Y(X, Y)
    X <- tidied$X
    Y <- tidied$Y
    if ((length(bg) > 0) & (is.vector(bg))) {
        bg <- matrix(bg, nrow = length(bg))
    }
    
    # select an epsilon (lowest non-zero value to use)
    epsilon <- min(Y[(Y > 0) & !is.na(Y)])
    
    
    # initial run to look for outliers:
    out0 <- deconLNR(
        Y = Y, X = X, bg = bg, weights = weights, epsilon = epsilon,
        maxit = maxit
    )
    # also get yhat and resids:
    out0$yhat <- X %*% out0$beta + bg
    out0$resids <- log2(pmax(Y, lower_thresh)) -
        log2(pmax(out0$yhat, lower_thresh))
    
    # ID bad genes:
    outliers <- flagOutliers(
        Y = Y,
        yhat = out0$yhat,
        wts = weights,
        resids = out0$resids,
        resid_thresh = resid_thresh
    )
    
    # remove outlier data points:
    Y.nooutliers <- replace(Y, outliers, NA)
    
    # re-run decon without outliers:
    out <- deconLNR(
        Y = Y.nooutliers,
        X = X,
        bg = bg,
        weights = weights,
        epsilon = epsilon
    )
    out$yhat <- X %*% out$beta + bg
    out$resids <- log2(pmax(Y.nooutliers, 0.5)) - log2(pmax(out$yhat, 0.5))
    
    # compute p-values
    tempbeta <- out$beta
    tempse <- tempp <- tempt <- tempbeta * NA
    for (i in seq_len(ncol(tempse))) {
        tempse[, i] <- suppressWarnings(sqrt(diag(out$sigmas[, , i])))
    }
    tempt <- (tempbeta / tempse)
    tempp <- 2 * (1 - stats::pnorm(tempt))
    out$p <- tempp
    out$t <- tempt
    out$se <- tempse
    
    # structure of output: beta, hessians, yhat, resids
    return(out)
}



#' Function to format Y, X inputs for decon
#'
#' Takes user-supplied X and Y, checks for accuracy, aligns by dimnames, adds
#' dimnames if missing
#'
#' @param X X matrix
#' @param Y Data matrix
#' @return X and Y, both formatted as matrices, with full dimnames and aligned
#' to each other by dimname
#' @keywords internal
#' @noRd
tidy_X_and_Y <- function(X, Y) {
    
    # format as matrices:
    Ynew <- Y
    if (is.vector(Y)) {
        Ynew <- matrix(Y, nrow = length(Y), dimnames = list(names(Y), "y"))
    }
    Xnew <- X
    
    # check alignment:
    if (!identical(rownames(Y), rownames(X))) {
        warning("Rows (genes) of X and Y are mis-aligned.")
    }
    out <- list(X = Xnew, Y = Ynew)
}
