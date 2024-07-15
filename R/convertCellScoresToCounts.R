
convertCellScoresToCounts <- function(beta, nuclei.counts = NULL,
                                      omit.tumor = FALSE) {
    # strip tumor rows if called for:
    if (omit.tumor) {
        beta <- beta[!grepl("tumor", rownames(beta)), , drop = FALSE]
    }
    
    # calc max abundance scores:
    max.total.abundance <- max(colSums(beta))
    
    # calculate rescaled scores:
    out <- list()
    out$cells.per.100 <- beta / max.total.abundance * 100
    if (length(nuclei.counts) == ncol(beta)) {
        out$cell.counts <- sweep(out$cells.per.100, 2, nuclei.counts, "*") / 100
    }
    return(out)
}
