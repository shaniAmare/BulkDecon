
download_profile_matrix <- function(species, age_group, matrixname) {
    
    # check formatting:
    if (length(species) > 1) {
        stop("specify just one species")
    }
    if (length(age_group) > 1) {
        stop("specify just one age")
    }
    if (length(matrixname) > 1) {
        stop("specify just one matrixname")
    }
    
    valid_species <- c("Human", "Mouse")
    
    if (!species %in% valid_species) {
        stop(paste0("Species input is invalid; must be \"", paste(valid_species, 
                                                                  collapse = "\" or \""), "\" (case sensitive)"))
    }
    
    if(species == "Human"){
        valid_ages <- c("Adult", "COVID-Infected", "Fetal")
    }else{
        valid_ages <- c("Adult", "Fetal/E14.5", "Fetal/E9.5-13.5", "Neonatal")
    }
    
    if (!age_group %in% valid_ages) {
        stop(paste0("Age input is invalid; must be \"", paste(valid_ages, collapse = "\" or \""), "\" (case sensitive)"))
    }
    
    metadata <- utils::read.csv(paste0("https://raw.github.com/Nanostring-Biostats/CellProfileLibrary/master/", species, "/",
                                       species, "_datasets_metadata.csv"), header = TRUE, sep = ",")
    
    librarynames <- paste0(metadata$Tissue, "_", metadata$Profile.Matrix)
    
    
    if (!is.element(matrixname, librarynames)) {
        stop(paste0(matrixname, " is not an expected cell profile matrix name for the species given. Did you mean \"", 
                    paste(librarynames[agrep(matrixname, librarynames)], collapse = "\" or \""), "\"?"))
    }
    
    profilename <- matrixname
    matrixname <- paste(species, age_group, matrixname, sep = "/")
    
    fulllibrarynames <- paste(metadata$Species, metadata$Age.Group, librarynames, sep = "/")
    
    if (!is.element(matrixname, fulllibrarynames)) {
        stop(paste0(profilename, " does not correspond with the age_group given. Did you mean \"", 
                    paste(unique(t(as.data.frame(strsplit(fulllibrarynames[agrep(profilename, fulllibrarynames)], "/")))[,2]), 
                          collapse = "\" or \""), "\"?"))
    }
    
    
    suppressMessages(repmis::source_data(paste0("https://raw.github.com/Nanostring-Biostats/CellProfileLibrary/master/", 
                                                matrixname, ".RData?raw=True"), 
                                         cache = FALSE, rdata = TRUE, envir = globalenv()))
    
    #assign("profile_matrix", as.matrix(profile_matrix), envir = globalenv())
    
    return(as.matrix(profile_matrix))
}

