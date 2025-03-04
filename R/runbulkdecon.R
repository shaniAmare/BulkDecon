
# idea here: 
# Use the bulk normallised count matrix - how?
# raw expression values (v$E from limma-voom)
# also use the single cell custom matrix that can be made
# check the anndata format into just using the mtrix
# run decon 

runbulkdecon <- function(X = NULL, norm_elt = NULL, raw_elt = NULL,
         wts = NULL,
         resid_thresh = 3, lower_thresh = 0.5,
         align_genes = TRUE,
         is_pure_tumor = NULL, n_tumor_clusters = 10,
         cell_counts = NULL,
         cellmerges = NULL,
         maxit = 1000){
  
  # check that norm and raw data exists:
  if(is.null(norm_elt)){
    stop("norm_elt must be set")
  }
  if(is.null(raw_elt)){
    stop("raw_elt must be set")
  }
  
  ##### HERE THE CHECKPOINT NEEDS TO BE TO IDENTIFY THE OBJECT TYPE FOR THE NORM AND RAW EXPRESSION DATA - probably not needed.
  # if (!is.element(norm_elt, names(object@assayData))) {
  #   stop(paste(norm_elt, "is not an element in assaysData slot"))
  # }
  # if (!is.element(raw_elt, names(object@assayData))) {
  #   stop(paste(raw_elt, "is not an element in assaysData slot"))
  # }
  
  # > object@assayData$q_norm %>% class()
  # [1] "matrix" "array" 
  # exprs is the same
  
  # prep components: - CHANGED TO ONLY HAVE only the MATRIX CONVERSION - checkpint or not?
  norm <- as.matrix(norm_elt)
  raw <- as.matrix( raw_elt)
  
  # estimate background
  bg <- calc_background(raw, norm)
  # bg <- derive_GeoMx_background(norm = norm,
  #                               # access the probe pool information from the feature metadata
  #                               probepool = Biobase::fData(object)$Module, # a character vector the length of the normalised matrix 
  #                               # access the names of the negative control probes 
  #                               ### TargetName IS GENE NAME, character vector and 
  #                               ### just extracting out of them what is called negative in the negative section of the fData.
  #                               # so negnames is also a character vector
  #                               negnames = Biobase::fData(object)$TargetName[Biobase::fData(object)$Negative])
  
  # run spatialdecon:
  result <- bulkdecon(norm = norm,
                         bg = bg,
                         X = X,
                         raw = raw,
                         wts = wts,
                         resid_thresh = resid_thresh, 
                         lower_thresh = lower_thresh,
                         align_genes = align_genes,
                         is_pure_tumor = is_pure_tumor, 
                         n_tumor_clusters = n_tumor_clusters,
                         cell_counts = cell_counts,
                         cellmerges = cellmerges,
                         maxit = maxit)
  return(result)
}