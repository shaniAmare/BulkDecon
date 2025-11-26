#' Collapse or Convert Cell Types Using a Linear Combination
#'
#' This function collapses individual cell-type rows of a matrix (typically
#' \code{beta} from a deconvolution fit) into higher-level groups using a
#' user-specified grouping structure. Optionally, covariance matrices
#' (\code{sigma}) are transformed accordingly.
#'
#' @param beta A numeric matrix with cell types in rows and samples in columns.
#' @param matching A named list. Each element is a character vector of
#'   row names from \code{beta} that should be grouped together. The list
#'   names define the collapsed output rows.
#' @param stat (Unused; reserved for future functionality.) Default: \code{sum}.
#' @param na.rm Logical, whether to remove NAs (unused; included for API
#'   compatibility).
#' @param sigma Optional covariance matrix or array. Can be:
#'   \itemize{
#'     \item a 2D covariance matrix (\eqn{k × k});
#'     \item a 3D array of covariances (\eqn{k × k × n}).
#'   }
#'
#' @details
#' The function constructs a linear transformation matrix \eqn{A} whose rows
#' correspond to collapsed cell-type groups and whose columns match the
#' original cell types. The collapsed abundance matrix is computed as:
#'
#' \deqn{beta' = A \, beta}
#'
#' If covariance matrices are supplied, they are transformed as:
#'
#' \deqn{Σ' = A \, Σ \, A^T}
#'
#' @return
#' If \code{sigma} is missing, returns only the collapsed \code{beta} matrix.
#'
#' If \code{sigma} is supplied, returns a list:
#' \describe{
#'   \item{beta}{Collapsed abundance matrix.}
#'   \item{sigma}{Collapsed covariance matrix or covariance array.}
#' }
#'
#' @export
convertCellTypes <- function(beta, matching, stat = sum,
                             na.rm = FALSE, sigma = NULL) {

  # ---------------- Input checks ------------------------------------------
  beta <- as.matrix(beta)

  if (!is.list(matching) || is.null(names(matching))) {
    stop("`matching` must be a *named* list of cell type groups.")
  }

  if (!all(unlist(matching) %in% rownames(beta))) {
    missing <- setdiff(unlist(matching), rownames(beta))
    stop("The following cell types in `matching` are missing in beta: ",
         paste(missing, collapse = ", "))
  }

  # ---------------- Build transformation matrix A -------------------------
  A <- matrix(
    0,
    nrow = length(matching),
    ncol = nrow(beta),
    dimnames = list(names(matching), rownames(beta))
  )

  for (grp in names(matching)) {
    A[grp, matching[[grp]]] <- 1
  }

  # ---------------- Apply transformation to beta --------------------------
  beta2 <- A %*% beta

  # ---------------- If no sigma, return beta2 only ------------------------
  if (is.null(sigma)) {
    return(beta2)
  }

  # ---------------- Transform sigma if provided ---------------------------
  if (length(dim(sigma)) == 2) {
    # Single covariance matrix
    sigma2 <- A %*% sigma %*% t(A)
  }
  else if (length(dim(sigma)) == 3) {
    # Covariance array per sample
    nsamp <- dim(sigma)[3]
    sigma2 <- array(
      NA,
      dim = c(nrow(A), nrow(A), nsamp),
      dimnames = list(rownames(A), rownames(A), dimnames(sigma)[[3]])
    )
    for (i in seq_len(nsamp)) {
      sigma2[, , i] <- A %*%
        sigma[, , i, drop = FALSE] %*%
        t(A)
    }
  }
  else {
    stop("`sigma` must be either a 2D covariance matrix or a 3D array.")
  }

  # ---------------- Return beta and sigma ---------------------------------
  return(list(beta = beta2, sigma = sigma2))
}
