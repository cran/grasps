#' Adjacency Matrix from Precision Matrix
#'
#' @description
#' Convert a precision matrix to a partial-correlation-based adjacency matrix.
#'
#' @param prec.mat A numeric precision matrix.
#'
#' @param diag.zero A logical value (default = TRUE) specifying whether to
#' set the diagonal entries of the adjacency matrix to 0.
#' If \code{diag.zero = FALSE}, the diagonal entries are set to 1 for a weighted
#' network. For an unweighted network (\code{weighted = FALSE}), the diagonal is
#' always forced to 0 to avoid self-loops.
#'
#' @param absolute A logical value (default = FALSE) specifying whether to
#' take the absolute values of the partial correlations.
#'
#' @param threshold A nonnegative numeric value (default = \code{NULL})
#' specifying the threshold for edge filtering.
#' Entries with absolute values smaller than the threshold are set to 0.
#'
#' @param weighted A logical value (default = TRUE) specifying whether to
#' return a weighted adjacency matrix.
#' If \code{weighted = FALSE}, the matrix is a binary adjacency matrix with
#' entries equal to 0 or 1.
#'
#' @return
#' A numeric adjacency matrix with S3 class \code{"adjmat"}.
#'
#' @details
#' For a precision matrix \eqn{\Omega}, the partial correlation between nodes
#' \eqn{i} and \eqn{j} is computed as
#' \deqn{\rho_{ij} = - \frac{\Omega_{ij}}{\sqrt{\Omega_{ii}\Omega_{jj}}}.}
#'
#' @example
#' inst/example/ex-grasps.R
#'
#' @export

prec_to_adj <- function(prec.mat, diag.zero = TRUE, absolute = FALSE,
                        threshold = NULL, weighted = TRUE) {

  if (!is.matrix(prec.mat)) {
    prec.mat <- as.matrix(prec.mat)
  }
  if (any(!is.finite(prec.mat))) {
    stop("The precision matrix contains non-finite values!")
  }
  if (!isSymmetric(prec.mat)) {
    stop("The precision matrix must be symmetric!")
  }
  d <- diag(prec.mat)
  if (any(d <= 0)) {
    stop("All diagonal entries of the precision matrix must be positive!")
  }

  ## partial-correlation-based adjacency matrix
  norm_mat <- outer(sqrt(d), sqrt(d), "*")
  adj_mat <- -prec.mat / norm_mat

  ## optional: diagonal entries
  if (diag.zero) {
    diag(adj_mat) <- 0
  } else {
    diag(adj_mat) <- 1
  }

  ## optional: absolute value
  if (absolute) {
    adj_mat <- abs(adj_mat)
  }

  ## optional: thresholding
  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1 || threshold < 0) {
      stop("The threshold must be a nonnegative scalar!")
    }
    adj_mat[abs(adj_mat) < threshold] <- 0
  }

  ## optional: binary adjacency matrix
  if (!weighted) {
    adj_mat <- 1 * (adj_mat != 0)
    diag(adj_mat) <- 0
    if (!diag.zero) {
      message("For an unweighted network, the diagonal is forced to 0 to avoid self-loops.")
    }
  }

  class(adj_mat) <- c("adjmat")
  return(adj_mat)
}
