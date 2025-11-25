bulkdecon <- function(norm,
                         bg,
                         X = NULL,
                         raw = NULL,
                         wts = NULL,
                         resid_thresh = 3, lower_thresh = 0.5,
                         align_genes = TRUE,
                         is_pure_tumor = NULL, n_tumor_clusters = 10,
                         cell_counts = NULL,
                         cellmerges = NULL,
                         maxit = 1000){

    #### preliminaries ---------------------------------

    # check formatting:
    if (!is.matrix(norm)) {
        stop("norm should be a matrix")
    }
    if ((length(X) > 0) & (!is.matrix(X))) {
        stop("X should be a matrix")
    }
    if ((length(raw) > 0) & (!is.matrix(raw))) {
        stop("raw must be a matrix")
    }
    if ((length(wts) > 0) & (!is.matrix(wts))) {
      stop("wts must be a matrix")
    }
    if ((length(cell_counts) > 0) & (!is.numeric(cell_counts))) {
        stop("cell_counts must be numeric")
    }


    if (length(bg) == 1) {
        bg <- matrix(bg, nrow(norm), ncol(norm),
                     dimnames = list(rownames(norm), colnames(norm))
        )
    }

    # If a matrix other than safeTME is input, rescale training matrix to avoid bad convergence properties:
    if (length(X) > 0) {
        # rescale matrix so its 99th percentile is near that of safeTME (which has a 99th percentile = 2.3)
        X <- X * 2 / quantile(X, 0.99)
    }

    # prep training matrix:
    if (length(X) == 0) {
        utils::data("safeTME", envir = environment())
        X <- BulkDecon::safeTME
    }
    sharedgenes <- intersect(rownames(norm), rownames(X))
    if (length(sharedgenes) == 0) {
        stop("no shared gene names between norm and X")
    }
    if (length(sharedgenes) < 100) {
        stop(paste0(
            "Only ", length(sharedgenes),
            " genes are shared between norm and X - this may not be enough
                to support accurate deconvolution."
        ))
    }

    # calculate weights based on expected SD of counts
    # wts = replace(norm, TRUE, 1)
    if (length(raw) > 0) {
        weight.by.TIL.resid.sd <-
            length(intersect(colnames(X), colnames(BulkDecon::safeTME))) > 10
        wts <- deriveWeights(norm,
                             raw = raw, error.model = "dsp",
                             weight.by.TIL.resid.sd = weight.by.TIL.resid.sd
        )
    }

    #### if pure tumor samples are specified, get tumor expression profile --------
    if (sum(is_pure_tumor) > 0) {

        # derive tumor profiles and merge into X:
        # (derive a separate profile for each tissue)
        X <- mergeTumorIntoX(
            norm = norm,
            bg = bg,
            pure_tumor_ids = is_pure_tumor,
            X = X[sharedgenes, ]
        )

        sharedgenes <- intersect(rownames(norm), rownames(X))
    }


    #### Run decon  -----------------------------------
    res <- algorithm2(
        Y = norm[sharedgenes, ],
        bg = bg[sharedgenes, ],
        X = X[sharedgenes, ],
        weights = wts[sharedgenes, ],
        maxit = maxit,
        resid_thresh = resid_thresh,
        lower_thresh = lower_thresh
    )


    #### combine closely-related cell types ------------------------------------

    if (length(cellmerges) > 0) {
        tempconv <- convertCellTypes(
            beta = res$beta,
            matching = cellmerges,
            stat = sum,
            na.rm = FALSE,
            sigma = res$sigmas
        )
        # overwrite original beta with merged beta:
        res$beta.granular <- res$beta
        res$sigma.granular <- res$sigmas
        res$sigmas <- NULL
        res$beta <- tempconv$beta
        res$sigma <- tempconv$sigma
    }


    #### compute p-values -------------------------------------------
    tempbeta <- res$beta
    tempse <- tempp <- tempt <- tempbeta * NA
    for (i in seq_len(ncol(tempse))) {
        tempse[, i] <- suppressWarnings(sqrt(diag(res$sigma[, , i])))
    }
    tempt <- (tempbeta / tempse)
    tempp <- 2 * (1 - stats::pnorm(tempt))
    res$p <- tempp
    res$t <- tempt
    res$se <- tempse


    #### rescale abundance estimates --------------------------------
    # (to proportions of total, proportions of immune, cell counts)

    # proportions:
    res$prop_of_all <- sweep(res$beta, 2, colSums(res$beta), "/")
    nontumorcellnames <- rownames(res$beta)[!grepl("tumor", rownames(res$beta))]
    res$prop_of_nontumor <- sweep(
        res$beta[nontumorcellnames, ], 2,
        colSums(res$beta[nontumorcellnames, ]), "/"
    )

    # on scale of cell counts:
    if (length(cell_counts) > 0) {
        res$cell.counts <- convertCellScoresToCounts(
            beta = res$beta,
            nuclei.counts = cell_counts,
            omit.tumor = TRUE
        )
        if (length(res$beta.granular) > 0) {
            res$cell.counts.granular <- convertCellScoresToCounts(
                beta = res$beta.granular,
                nuclei.counts = cell_counts,
                omit.tumor = TRUE
            )
        }
    }

    # add other pertinent info to res:
    res$X <- X[rownames(res$resids), ]
    return(res)
}
