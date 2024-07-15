
mergeTumorIntoX <- function(norm, bg, pure_tumor_ids, X) {
    
    # round up 0 values in norm:
    min.nonzero <- min(norm[norm > 0], na.rm = TRUE)
    norm <- pmax(norm, min.nonzero)
    
    # subset data to only the pure tumor IDs:
    norm <- norm[, pure_tumor_ids, drop = FALSE]
    bg <- bg[, pure_tumor_ids, drop = FALSE]
    
    # bg-subtract:
    norm <- pmax(norm - bg, min(norm) / 20)
    
    # # fix K if too big:
    # if (ncol(norm) < K) {
    #     K <- ncol(norm)
    # }
    # 
    # # case 1: want to use every column in norm as a separate profile:
    # #  (includes case of just one column in norm)
    # if (K == ncol(norm)) {
    #     tumorX <- norm
    # }
    # 
    # # case 2: if many tumor AOIs, get profiles for K clusters of data:
    # if (K < ncol(norm)) {
    #     # cluster and cut:
    #     h <- stats::hclust(stats::dist(t(log2(norm))))
    #     cut <- stats::cutree(h, k = K)
    #     # get clusters' geomean profiles:
    #     tumorX <- c()
    #     for (cid in unique(cut)) {
    #         tumorX <- cbind(
    #             tumorX,
    #             exp(rowMeans(log(norm[, cut == cid, drop = FALSE])))
    #         )
    #     }
    #     colnames(tumorX) <- paste0("tumor.", seq_len(ncol(tumorX)))
    # }
    # 
    # align norm with X:
    sharedgenes <- intersect(rownames(norm), rownames(X))
    tumorX <- norm[sharedgenes, , drop = FALSE]
    X <- X[sharedgenes, , drop = FALSE]
    
    # rescale tumor X:
    meanq90 <- max(mean(apply(X, 2, stats::quantile, 0.9)), 1e-3)
    tumorq90s <- apply(tumorX, 2, stats::quantile, 0.9)
    tumorX <- sweep(tumorX, 2, tumorq90s, "/") * meanq90
    
    # merge:
    out <- cbind(X, tumorX)
    return(out)
}
