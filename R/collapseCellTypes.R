#' Collapse Cell Types into Higher-Level Groups
#'
#' This function collapses cell-type–level deconvolution results into
#' broader categories by summing selected rows of the coefficient matrix
#' (`beta`) and—if present—related quantities such as
#' `prop_of_all`, `prop_of_nontumor`, and covariance matrices (`sigmas`).
#'
#' @param fit A deconvolution results list, typically containing
#'   \code{beta}, and optionally \code{prop_of_all}, \code{prop_of_nontumor},
#'   \code{sigmas}, and \code{X}.
#' @param matching A named list where each element is a character vector of
#'   cell-type names to be collapsed. The names of the list become the
#'   final collapsed categories.
#'
#' @details
#' The function creates a linear transformation matrix \eqn{A} that
#' aggregates the specified cell types. For example:
#' \preformatted{
#' matching = list(
#'     Immune = c("T_cells", "B_cells"),
#'     Stromal = c("Fibroblasts", "Endothelium")
#' )
#' }
#'
#' produces a 2 × k matrix \eqn{A}, which is applied as:
#' \deqn{beta_collapsed = A \times beta}
#'
#' If covariance arrays (`sigmas`) are supplied (2D or 3D), they are
#' transformed as:
#' \deqn{Σ' = A Σ A^T}
#'
#' @return
#' A modified copy of \code{fit} in which the relevant matrices
#' have been collapsed. Returns:
#' \describe{
#'   \item{beta}{Collapsed abundance matrix.}
#'   \item{prop_of_all}{Collapsed if present.}
#'   \item{prop_of_nontumor}{Collapsed if present.}
#'   \item{sigmas}{Collapsed covariance matrices if present.}
#'   \item{se, t, p}{Recomputed standard errors, t-stats, and p-values
#'        if covariance matrices are available.}
#' }
#'
#' @export
collapseCellTypes <- function(fit, matching) {

  # ---- Input checks -------------------------------------------------------
  if (is.null(fit) || !is.list(fit)) {
    stop("`fit` must be a list produced by a deconvolution function.")
  }
  if (!is.list(matching) || is.null(names(matching))) {
    stop("`matching` must be a *named* list of cell-type groups.")
  }
  if (!"beta" %in% names(fit)) {
    stop("`fit` must contain a `beta` matrix.")
  }

  beta <- fit$beta
  beta <- as.matrix(beta)

  # Ensure all referenced cell types exist
  all_cells <- rownames(beta)
  requested_cells <- unique(unlist(matching))
  missing <- setdiff(requested_cells, all_cells)

  if (length(missing) > 0) {
    stop("The following cell types in `matching` are not in beta: ",
         paste(missing, collapse = ", "))
  }

  # ---- Build linear transformation matrix A -------------------------------
  startingcellnames <- requested_cells
  A <- matrix(
    0,
    nrow = length(matching),
    ncol = length(all_cells),
    dimnames = list(names(matching), all_cells)
  )

  for (grp in names(matching)) {
    A[grp, matching[[grp]]] <- 1
  }

  # ---- Begin output -------------------------------------------------------
  out <- fit

  # ---- Apply A to beta and other abundance matrices -----------------------
  to_collapse <- c("beta", "prop_of_all", "prop_of_nontumor")

  for (nm in to_collapse) {
    if (nm %in% names(fit)) {
      mat <- as.matrix(fit[[nm]])
      out[[nm]] <- A[, startingcellnames, drop = FALSE] %*%
        mat[startingcellnames, , drop = FALSE]
    }
  }

  # ---- Collapse covariance matrices if present ----------------------------
  if ("sigmas" %in% names(fit)) {
    sigma <- fit$sigmas

    # Case 1: 2D matrix (single covariance)
    if (length(dim(sigma)) == 2) {
      sigma_sub <- sigma[startingcellnames, startingcellnames, drop = FALSE]
      out$sigmas <- A[, startingcellnames, drop = FALSE] %*%
        sigma_sub %*% t(A[, startingcellnames, drop = FALSE])
    }

    # Case 2: 3D array (per-sample covariance)
    if (length(dim(sigma)) == 3) {
      nsamp <- dim(sigma)[3]
      A_sub <- A[, startingcellnames, drop = FALSE]

      out$sigmas <- array(
        NA,
        dim = c(nrow(A), nrow(A), nsamp),
        dimnames = list(rownames(A), rownames(A), dimnames(sigma)[[3]])
      )

      for (i in seq_len(nsamp)) {
        sig_i <- sigma[startingcellnames, startingcellnames, i, drop = FALSE]
        out$sigmas[, , i] <- A_sub %*% sig_i %*% t(A_sub)
      }
    }
  }

  # ---- Recompute SE, t-stats, p-values -----------------------------------
  if ("beta" %in% names(out) && "sigmas" %in% names(out)) {

    tempbeta <- out$beta
    n_groups <- nrow(tempbeta)
    n_samples <- ncol(tempbeta)

    tempse <- matrix(NA, nrow = n_groups, ncol = n_samples,
                     dimnames = dimnames(tempbeta))

    # Extract SE per-sample
    if (length(dim(out$sigmas)) == 2) {
      diag_se <- sqrt(diag(out$sigmas))
      tempse[,] <- diag_se
    } else {
      for (i in seq_len(n_samples)) {
        tempse[, i] <- suppressWarnings(
          sqrt(diag(out$sigmas[, , i]))
        )
      }
    }

    out$se <- tempse
    out$t <- out$beta / tempse

    # Only compute p-values if original design matrix existed
    if ("X" %in% names(out)) {
      out$p <- 2 * (1 - stats::pnorm(abs(out$t)))
    }
  }

  return(out)
}
