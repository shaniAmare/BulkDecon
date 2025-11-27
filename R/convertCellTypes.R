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
  # format matching list as a matrix to take a linear combination of beta:
  A <- matrix(0, length(matching), nrow(beta),
              dimnames = list(names(matching), rownames(beta))
  )
  for (name in names(matching)) {
    cellnames <- matching[[name]]
    A[name, cellnames] <- 1
  }

  # apply A transformation to beta:
  beta2 <- A %*% beta

  # if Sigma provided, get vcov of beta2:
  if (length(sigma) > 0) {
    if (length(dim(sigma)) == 2) {
      sigma2 <- A %*% sigma %*% t(A)
    }
    if (length(dim(sigma)) == 3) {
      sigma2 <- array(NA,
                      dim = c(nrow(A), nrow(A), dim(sigma)[3]),
                      dimnames = list(rownames(A), rownames(A), dimnames(sigma)[[3]])
      )
      for (i in seq_len(dim(sigma)[3])) {
        sigma2[, , i] <- A %*% sigma[, , i] %*% t(A)
      }
    }
  }

  # if no Sigma, just return transformed beta:
  if (length(sigma) == 0) {
    return(beta2)
  }
  # if there is a sigma, return beta and the sigma:
  if (length(sigma) > 0) {
    out <- list(beta = beta2, sigma = sigma2)
    return(out)
  }
}
