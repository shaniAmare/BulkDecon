deriveWeights <- function(norm, raw = NULL, error.model = "dsp",
                          weight.by.TIL.resid.sd = FALSE) {

    # get tech SDs if raw data provided:
    if (length(raw) == 0) {
        sds.tech <- matrix(0.1, nrow(raw), ncol(raw), dimnames = dimnames(raw))
    }
    if (length(raw) > 0) {
        sds.tech <- runErrorModel(
            counts = raw,
            platform = "dsp"
        )
    }

    # if the mean.resid.sd vector (which defines genes' biological SD) is in
    # the environment, get biological noise:
    if (!weight.by.TIL.resid.sd) {
        sds.bio <- matrix(0.1, nrow(raw), ncol(raw), dimnames = dimnames(raw))
    }
    if (weight.by.TIL.resid.sd) {
        utils::data("mean.resid.sd", envir = environment())
        sds.bio <- matrix(NA, nrow(raw), ncol(raw), dimnames = dimnames(raw))
        for (gene in intersect(
            names(BulkDecon::mean.resid.sd),
            rownames(sds.bio)
        )) {
            sds.bio[gene, ] <- BulkDecon::mean.resid.sd[gene]
        }
        sds.bio <- replace(sds.bio, is.na(sds.bio), mean(sds.bio, na.rm = TRUE))
    }

    # define total SD, and invert to get weights
    sds.tot <- sqrt(sds.tech^2 + sds.bio^2)
    wts <- 1 / sds.tech
    return(wts)
}
