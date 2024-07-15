
collapseCellTypes <- function(fit, matching) {
    
    # results object to hold the collapsed results:
    out <- fit
    
    # format matching list as a matrix to take a linear combination of beta:
    startingcellnames <- unlist(matching)
    A <- matrix(0, length(matching), nrow(fit$beta),
                dimnames = list(names(matching), rownames(fit$beta))
    )
    for (name in names(matching)) {
        cellnames <- matching[[name]]
        A[name, cellnames] <- 1
    }
    
    # apply A transformation to beta:
    for (name in c("beta", "prop_of_all", "prop_of_nontumor")) {
        if (is.element(name, names(fit))) {
            out[[name]] <- A[, startingcellnames] %*% fit[[name]][startingcellnames, ]
        }
    }
    
    # if Sigma provided, get vcov of beta2:
    if (is.element("sigmas", names(out))) {
        sigma <- fit$sigmas
        if (length(dim(sigma)) == 2) {
            out$sigmas <- A[, startingcellnames] %*%
                sigma[startingcellnames, startingcellnames, ] %*%
                t(A[, startingcellnames])
        }
        if (length(dim(sigma)) == 3) {
            out$sigmas <- array(NA,
                                dim = c(nrow(A), nrow(A), dim(sigma)[3]),
                                dimnames = list(rownames(A), rownames(A), dimnames(sigma)[[3]])
            )
            for (i in seq_len(dim(sigma)[3])) {
                out$sigmas[, , i] <- A %*% sigma[, , i] %*% t(A)
            }
        }
    }
    
    # re-calculate p, se, t:
    if (is.element("beta", names(out)) & is.element("sigmas", names(out))) {
        # compute p-values
        tempbeta <- out$beta
        tempse <- tempp <- tempt <- tempbeta * NA
        for (i in seq_len(ncol(tempse))) {
            tempse[, i] <- suppressWarnings(sqrt(diag(out$sigmas[, , i])))
        }
        out$se <- tempse
        out$t <- (tempbeta / tempse)
        if (is.element("X", names(out))) {
            out$p <- 2 * (1 - stats::pnorm(out$t))
        }
    }
    
    return(out)
}
