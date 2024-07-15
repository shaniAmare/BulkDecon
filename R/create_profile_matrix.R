create_profile_matrix <- function(mtx, cellAnnots, cellTypeCol, cellNameCol,
                                  matrixName = "Custom", outDir = "./", 
                                  geneList = NULL, normalize = FALSE,
                                  scalingFactor = 5, minCellNum = 15,
                                  minGenes = 100, discardCellTypes = FALSE) {
    
    # checking user input values
    if(is.null(mtx)){
        stop("count matrix is necessary")
    }
    if(!is.null(outDir)){
        if (!dir.exists(outDir) ) {
            stop("Output directory is not valid")
        }
    }
    if (is.null(cellAnnots)){
        stop("Cell Annotations are needed")
    }
    
    if(is.null(cellTypeCol)){
        stop("cellTypeCol must not be NULL")
    }else if (!cellTypeCol %in% colnames(cellAnnots)){
        stop("cellTypeCol not in cellAnnots")
    }
    
    if(is.null(cellNameCol)){
        stop("cellNameCol must not be NULL")
    }else if (!cellNameCol %in% colnames(cellAnnots)){
        stop("cellNameCol not in cellAnnots")
    }
    
    if(!is(normalize, "logical")){
        warning("normalize not a boolean, continuing with assumption that data should not be normalized")
        normalize <- TRUE
    }
    if(!is(discardCellTypes, "logical")){
        warning("discardCellTypes not a boolean, continuing with default of discarding cell types")
        discardCellTypes <- TRUE
    }
    
    if(!is(scalingFactor, "numeric")){
        warning("scalingFactor not a numeric, continuing with default value of 5")
        scalingFactor <- 5
    }
    if(!is(minCellNum, "numeric")){
        warning("minCellNum not a numeric, continuing with default value of 15")
        minCellNum <- 15
    }
    if(!is(minGenes, "numeric")){
        warning("minGenes not a numeric, continuing with default value of 100")
        minGenes <- 100
    }
    
    #make a sparse matrix
    mtx <- Matrix::Matrix(as.matrix(mtx), sparse = TRUE) 
    
    cellTypes <- NULL
    
    #read in cell type annotation file
    #get cell types 
    cellTypes <- cellAnnots[[cellTypeCol]]
    #assign cell name to type
    names(cellTypes) <- cellAnnots[[cellNameCol]]
    
    if(is.null(cellTypes)){
        stop("cellAnnots and/or cellTypeCol arguments are incorrectly formatted. 
                 cellAnnots should be a data frame, and cellTypeCol should give the 
                 name of the column holding each cell's cell type.")
    }
    
    rm(cellAnnots)
    
    # mtx <- as.data.frame(mtx)
    
    if(!any(names(cellTypes) %in% colnames(mtx)) & 
       any(names(cellTypes) %in% rownames(mtx))){
        print("Transposing Matrix")
        mtx<- t(mtx)
    }
    
    if(!any(names(cellTypes) %in% colnames(mtx))){
        stop(paste("cellNameCol names does not match count matrix column names", 
                   "matrix cell names:", colnames(mtx)[1], "annots cell names:", names(cellTypes)[1]))
    }else if(!all(names(cellTypes) %in% colnames(mtx))){
        missing <- length(which(!names(cellTypes) %in% colnames(mtx)))
        warning(paste("not all cellNameCol names are in count matrix;", missing, "cells are missing"))
    }
    
    if(discardCellTypes == TRUE){
        #remove cells with no cell type assignment
        w2rm <- which(is.na(cellTypes) | tolower(cellTypes) %in% c("unspecified", "unknown", "not available")) 
        w2rm <- unique(c(w2rm, grep(pattern = "doublet|dividing|low q|filtered|mitotic", x = tolower(cellTypes))))
        if(length(w2rm) > 0){
            cellTypes <- cellTypes[-w2rm]
        }
    }
    
    #normalize data if necessary 
    if(normalize == TRUE){
        print("Normalizing Matrix")
        med <- median(Matrix::colSums(mtx))
        cols <- colnames(mtx)
        rows <- rownames(mtx)
        mtx <- Matrix::Matrix(sweep(mtx, 2, Matrix::colSums(mtx), "/") * med, 
                              sparse = TRUE) 
        
        colnames(mtx) <- cols
        rownames(mtx) <- rows
        
        rm(cols,rows)
    }
    
    atlas <- NULL
    
    #get all unique cell types
    CTs <- unique(cellTypes)
    
    print("Creating Atlas")
    
    #change to apply() if bioconductor requires it
    for(i in CTs){
        #print log of progress
        print(paste(which(CTs == i), "/", length(CTs), ":", i))
        
        #get cell names for this cell type
        cellsType <- names(cellTypes)[which(cellTypes == i)]
        #confirm cells are in matrix
        cellsType <- cellsType[which(cellsType %in% colnames(mtx))]
        
        if(length(cellsType) > minCellNum){
            
            if(length(cellsType) > 1){
                #remove cells with low gene expression
                cellsType <- cellsType[which(Matrix::colSums(mtx[,cellsType] > 0) > minGenes)]
            }else{
                cellsType <- cellsType[which(sum(mtx[,cellsType] > 0) > minGenes)]
            }
            
            if(length(cellsType) > minCellNum){
                #get average expression if there are enough cells for cell type
                if(length(cellsType) > minCellNum & length(cellsType) != 1){
                    atlas <- as.data.frame(cbind(atlas, Matrix::rowMeans(mtx[,cellsType], na.rm = TRUE)))
                    colnames(atlas)[ncol(atlas)] <- i
                }else{
                    atlas <- as.data.frame(cbind(atlas, mtx[,cellsType]))
                    colnames(atlas)[ncol(atlas)] <- i
                }
            }else{
                warning(paste("\n", i, "was dropped from matrix because it didn't have enough viable cells based on current filtering thresholds. 
                                        If this cell type is necessary consider changing minCellNum or minGenes\n"))
            }
        }else{
            warning(paste("\n", i, "was dropped from matrix because it didn't have enough viable cells based on current filtering thresholds. 
                                        If this cell type is necessary consider changing minCellNum or minGenes\n"))
        }
    }
    
    rm(cellTypes)
    
    numCellTypesExpr <- 1
    
    #subset to genes expressed in at least a user defined number of cell type(s)
    if(ncol(atlas) == 1){
        w2kp <- which(Matrix::rowSums(atlas > 0) >= numCellTypesExpr)
        cols <- colnames(atlas)
        rows <- rownames(atlas)[w2kp]
        
        atlas <- as.matrix(atlas[w2kp,])
        
        colnames(atlas) <- cols
        rownames(atlas) <- rows
        
        rm(cols, rows, w2kp)
    }else{
        atlas <- atlas[which(Matrix::rowSums(atlas > 0) >= numCellTypesExpr),]
    }
    
    #scale data
    atlas <- atlas * scalingFactor
    
    if(!is.null(geneList)){
        if(any(geneList %in% rownames(atlas))){
            #subset to genes in panel
            atlas <- atlas[rownames(atlas) %in% geneList,]
        }else{
            warning("geneList genes do not match genes in matrix, no filtering done")
        }
    }
    
    #ensure cell types don't contain ","
    colnames(atlas) <- gsub(pattern = ",", replacement = "-", x = colnames(atlas))
    
    if(!is.null(outDir)){
        #write profile matrix
        write.table(atlas, file = paste0(outDir, "/", matrixName, "_profileMatrix.csv"), 
                    row.names = TRUE, col.names = NA, quote = FALSE, sep = ",")
    }
    
    return(as.matrix(atlas))
}

