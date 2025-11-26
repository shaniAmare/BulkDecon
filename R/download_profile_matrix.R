#' Download a Cell Profile Matrix from the NanoString CellProfileLibrary
#'
#' @description
#' Downloads a cell profile matrix (RData format) from the public NanoString
#' CellProfileLibrary GitHub repository, based on `species`, `age_group`,
#' and `matrixname`.
#'
#' @param species Character string. One of `"Human"` or `"Mouse"`.
#' @param age_group Character string. Valid options depend on the species:
#'   \itemize{
#'     \item Human: `"Adult"`, `"COVID-Infected"`, `"Fetal"`
#'     \item Mouse: `"Adult"`, `"Fetal/E14.5"`, `"Fetal/E9.5-13.5"`, `"Neonatal"`
#'   }
#' @param matrixname Character string giving the desired profile matrix name
#'   (e.g. `"Liver_Immune"`, `"Colon_Fibroblast"`).
#'
#' @return A numeric matrix containing the downloaded cell profile matrix.
#'
#' @export
download_profile_matrix <- function(species, age_group, matrixname) {

  ## 1. Basic input checks ------------------------------------------------------
  if (length(species) != 1) stop("Specify exactly one species.")
  if (length(age_group) != 1) stop("Specify exactly one age_group.")
  if (length(matrixname) != 1) stop("Specify exactly one matrixname.")

  valid_species <- c("Human", "Mouse")

  if (!species %in% valid_species) {
    stop(paste0("Invalid species. Must be one of: ",
                paste(valid_species, collapse = ", "),
                " (case sensitive)."))
  }

  valid_ages <- if (species == "Human") {
    c("Adult", "COVID-Infected", "Fetal")
  } else {
    c("Adult", "Fetal/E14.5", "Fetal/E9.5-13.5", "Neonatal")
  }

  if (!age_group %in% valid_ages) {
    stop(paste0("Invalid age_group for species ", species,
                ". Must be one of: ",
                paste(valid_ages, collapse = ", "),
                " (case sensitive)."))
  }

  ## 2. Load metadata -----------------------------------------------------------
  metadata_url <- paste0(
    "https://raw.githubusercontent.com/Nanostring-Biostats/CellProfileLibrary/master/",
    species, "/", species, "_datasets_metadata.csv"
  )

  metadata <- utils::read.csv(metadata_url, stringsAsFactors = FALSE)

  librarynames <- paste0(metadata$Tissue, "_", metadata$Profile.Matrix)

  ## 3. Validate matrixname -----------------------------------------------------
  if (!matrixname %in% librarynames) {
    suggestion <- librarynames[agrep(matrixname, librarynames)]
    stop(paste0(
      matrixname, " is not a valid matrix name for ", species, ". ",
      if (length(suggestion))
        paste0("Did you mean: ", paste(suggestion, collapse = " or "), "?")
      else
        ""
    ))
  }

  ## 4. Validate species-age-matrix combo --------------------------------------
  fullnames <- paste(metadata$Species, metadata$Age.Group, librarynames, sep = "/")

  requested <- paste(species, age_group, matrixname, sep = "/")

  if (!requested %in% fullnames) {
    # suggest available ages for this matrix
    match_idx <- agrep(matrixname, fullnames)
    suggested_ages <- unique(
      vapply(strsplit(fullnames[match_idx], "/"), `[`, "", 2)
    )

    stop(paste0(
      "Selected matrix '", matrixname,
      "' does not match age_group '", age_group, "'. ",
      "Available age groups for this matrix: ",
      paste(suggested_ages, collapse = ", "), "."
    ))
  }

  ## 5. Build download URL ------------------------------------------------------
  file_url <- paste0(
    "https://raw.githubusercontent.com/Nanostring-Biostats/CellProfileLibrary/master/",
    requested, ".RData?raw=True"
  )

  ## 6. Load RData object -------------------------------------------------------
  tmpfile <- tempfile(fileext = ".RData")
  utils::download.file(file_url, destfile = tmpfile, quiet = TRUE)

  e <- new.env()
  load(tmpfile, envir = e)

  if (!"profile_matrix" %in% ls(e)) {
    stop("Downloaded RData file does not contain an object named 'profile_matrix'.")
  }

  return(as.matrix(e$profile_matrix))
}
