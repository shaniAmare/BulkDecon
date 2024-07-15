
flagOutliers <- function(Y, yhat, resids, wts, resid_thresh = 3) {
    
    # get weighted resids:
    if (length(wts) == 0) {
        wres <- resids
    }
    if (length(wts) > 0) {
        wres <- resids * wts
    }
    
    # flag bad genes:
    outlier_genes <- c() # <-- this line makes it so no outlier genes are filtered
    # flag bad data points: (not doing anything for now)
    outlier_data_points <- abs(resids) > resid_thresh
    return(outlier_data_points)
}
