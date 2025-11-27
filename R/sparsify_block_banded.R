#' Groupwise Block-Banded Sparsifier
#'
#' @description
#' Make a precision-like matrix block-banded according to group membership,
#' keeping only entries within specified group neighborhoods.
#'
#' @param mat A \eqn{d \times d} precision-like matrix specifying the base
#' matrix to be masked.
#'
#' @param membership An integer vector specifying the group membership.
#' The length of \code{membership} must be consistent with the dimension \eqn{d}.
#'
#' @param neighbor.range An integer (default = 1) specifying the neighbor range,
#' where groups whose labels differ by at most \code{neighbor.range} are
#' considered neighbors and kept in the mask.
#'
#' @return An object with S3 class \code{"sparsify_block_banded"} containing
#' the following components:
#' \describe{
#' \item{Omega}{The masked precision matrix.}
#' \item{Sigma}{The covariance matrix, i.e., the inverse of \code{Omega}.}
#' \item{sparsity}{Proportion of zero entries in \code{Omega}.}
#' \item{membership}{An integer vector specifying the group membership.}
#' }
#'
#' @example
#' inst/example/ex-sparsify_block_banded.R
#'
#' @export

sparsify_block_banded <- function(mat, membership, neighbor.range = 1) {

  d <- ncol(mat)
  if (length(membership) != d) {
    stop(sprintf("Length of 'membership' (%d) must equal the matrix dimension d (%d).",
                 length(membership), d))
  }

  ## determine which entries to keep: membership within 'neighbor.range'
  mask <- abs(outer(membership, membership, `-`)) <= neighbor.range

  ## apply mask to matrix
  Omega <- mat * mask
  ## compute covariance
  Sigma <- solve(Omega)

  result <- list(Omega = Omega, Sigma = Sigma,
                 sparsity = sum(Omega == 0) / length(Omega),
                 membership = membership)
  class(result) <- c("sparsify_block_banded", "blkmat")
  return(result)
}
