
calc_background <- function(exprs_mat, norm){ #, norm, probepool, negnames) {
    
    # # check data input:
    # if (nrow(norm) != length(probepool)) {
    #     stop("nrow(norm) != length(probepool)")
    # }
    
    # initialize:
    bg <- exprs_mat * 0
    
    results <- calc_bgindex(exprs_mat)
    # for above function (i.e., calc_bgindex)
    # Genes with a stability index rank percentile >80 as well as a reversed rank percentile >60 for each of the 4 stability feature
    # The main statistic is the segIdx column, which is the SEG index - and I changed it to bgIdx
    negnames <- dplyr::filter_at(results, vars(bgIdx), ~. > quantile(., probs = 0.8)) %>% row.names()
    # fill in expected background at scale of normalized data:
    # for (pool in unique(probepool)) {
    #     
    tempnegs <- intersect(negnames, rownames(exprs_mat))
    if (length(tempnegs) == 0) {
        stop(paste0(pool, " did not identify any background stable genes. did you input the raw bulk count matrix?"))
    }
    # following is a numeric vector()
    tempnegfactor <- colMeans(norm[tempnegs, , drop = FALSE])
    
    # fill in the corresponding elements of bg:
    bg<- sweep(bg, 2, tempnegfactor, "+")
    # }
    return(bg)
}
