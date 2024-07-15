
convertCellTypes <- function(beta, matching, stat = sum,
                             na.rm = FALSE, sigma = NULL) {
    # format matching list as a matrix to take a linear combination of beta:
    A <- matrix(0, length(matching), nrow(beta),
                dimnames = list(names(matching), rownames(beta))
    )
    for (name in names(matching)) {
        cellnames <- matching[[name]]
        A[name, cellnames] <- 1
    }
    
    # apply A transformation to beta:
    beta2 <- A %*% beta
    
    # if Sigma provided, get vcov of beta2:
    if (length(sigma) > 0) {
        if (length(dim(sigma)) == 2) {
            sigma2 <- A %*% sigma %*% t(A)
        }
        if (length(dim(sigma)) == 3) {
            sigma2 <- array(NA,
                            dim = c(nrow(A), nrow(A), dim(sigma)[3]),
                            dimnames = list(rownames(A), rownames(A), dimnames(sigma)[[3]])
            )
            for (i in seq_len(dim(sigma)[3])) {
                sigma2[, , i] <- A %*% sigma[, , i] %*% t(A)
            }
        }
    }
    
    # if no Sigma, just return transformed beta:
    if (length(sigma) == 0) {
        return(beta2)
    }
    # if there is a sigma, return beta and the sigma:
    if (length(sigma) > 0) {
        out <- list(beta = beta2, sigma = sigma2)
        return(out)
    }
}
