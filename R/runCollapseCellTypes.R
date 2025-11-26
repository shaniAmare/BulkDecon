#' Collapse cell types in a deconvolution result object
#'
#' Wrapper that extracts deconvolution outputs from a Biobase::ExpressionSet,
#' applies \code{collapseCellTypes()}, and writes the collapsed results back
#' into the object's \code{pData}.
#'
#' @param object An \code{ExpressionSet} containing deconvolution outputs
#'   stored in its \code{pData} slot (beta, p, t, se, prop_of_all,
#'   prop_of_nontumor, sigmas).
#' @param matching A named list specifying how fine-grained cell types should
#'   be collapsed into broader groups. Required.
#'
#' @return The input \code{ExpressionSet} with collapsed cell-type results
#'   written back into \code{pData}.
#' @export
runCollapseCellTypes <- function(object, matching = NULL) {

  # --- Input checks ---------------------------------------------------------
  if (is.null(matching)) {
    stop("`matching` must be provided and cannot be NULL.")
  }

  if (!inherits(object, "ExpressionSet")) {
    stop("`object` must be a Biobase::ExpressionSet.")
  }

  pdat <- Biobase::pData(object)

  required_fields <- c(
    "beta", "p", "t", "se",
    "prop_of_all", "prop_of_nontumor", "sigmas"
  )

  missing_fields <- setdiff(required_fields, colnames(pdat))
  if (length(missing_fields) > 0) {
    stop("The following required fields are missing from pData(object): ",
         paste(missing_fields, collapse = ", "))
  }

  # --- Build input list for collapseCellTypes -------------------------------
  fit <- list(
    beta            = t(pdat$beta),
    p               = t(pdat$p),
    t               = t(pdat$t),
    se              = t(pdat$se),
    prop_of_all     = t(pdat$prop_of_all),
    prop_of_nontumor = t(pdat$prop_of_nontumor),
    sigma           = pdat$sigmas
  )

  # --- Run collapse ---------------------------------------------------------
  temp <- collapseCellTypes(fit = fit, matching = matching)

  # --- Write collapsed results back ----------------------------------------
  Biobase::pData(object)$beta            <- t(temp$beta)
  Biobase::pData(object)$p               <- t(temp$p)
  Biobase::pData(object)$t               <- t(temp$t)
  Biobase::pData(object)$se              <- t(temp$se)
  Biobase::pData(object)$prop_of_all     <- t(temp$prop_of_all)
  Biobase::pData(object)$prop_of_nontumor <- t(temp$prop_of_nontumor)
  Biobase::pData(object)$sigmas          <- temp$sigma

  return(object)
}
