#' Line search for the sharp upper bound lambda_max
#'
#' @description
#' Given a list of off-diagonal blocks \code{blockList} from the sample
#' covariance (or a related matrix), this routine finds a \strong{sharp}
#' \eqn{\lambda_{\max}} by bisection so that the KKT condition
#' \deqn{\Vert\text{soft}(B, \lambda\alpha)\Vert_F \leq \lambda(1-\alpha)}
#' holds for every block \eqn{B} in \code{blockList}, where
#' \eqn{\text{soft}(\cdot, t)} denotes element-wise soft-thresholding with
#' threshold \eqn{t}.
#'
#' @param blockList A list of numeric matrices where each element is a block
#' used in the group KKT check.
#'
#' @param lambda.safe A numeric value representing a safe (but possibly loose)
#' upper bound for \eqn{\lambda_{\max}}. If the value is negative, the search
#' starts from 1.
#'
#' @param alpha A numeric value in (0, 1) specifying the mixture weight between
#' the individual (L1) and group (L2) penalties in sparse-group formulations.
#'
#' @param growiter An integer specifying the maximum number of exponential
#' growth steps during the initial search for an admissible upper bound
#' \eqn{\lambda_{\max}}.
#'
#' @param tol A numeric value > 0 specifying the relative tolerance for
#' the bisection stopping rule on the interval width.
#'
#' @param maxiter An integer specifying the maximum number of bisection
#' iterations.
#'
#' @return
#' A numeric value giving the bisection-refined upper bound \eqn{\lambda_{\max}}.
#'
#' @details
#' The procedure first grows an upper bound starting from \code{lambda.safe}
#' (or 1 if that is non-positive), doubling until the KKT check passes or
#' a growth cap is reached, and then performs bisection between the current
#' lower/upper bounds until the interval width is below
#' \code{tol * max(1, lambda_upper)} or \code{maxiter} is hit.
#'
#' @noRd

line_search_lambda_max <- function(blockList, lambda.safe, alpha, growiter, tol, maxiter) {
  stopifnot(is.list(blockList), length(blockList) > 0, alpha > 0, alpha < 1)
  KKTcheck <- function(lambda) {
    for (B in blockList) {
      Thres <- matrix(lambda*alpha, nrow(B), ncol(B))
      if (norm(soft_matrix(B, Thres), "F") > lambda*(1-alpha)+tol) {
        return(FALSE)
      }
    }
    return(TRUE)
  }
  lambda_upper <- if (lambda.safe > 0) lambda.safe else 1
  grow <- 0
  while (!KKTcheck(lambda_upper) && grow < growiter) {
    lambda_upper <- 2 * lambda_upper
    grow <- grow + 1
  }
  lambda_lower <- 0
  iter <- 0
  while (iter < maxiter && (lambda_upper-lambda_lower) > tol*max(1,lambda_upper)) {
    iter <- iter + 1
    lambda_mid <- (lambda_lower + lambda_upper) / 2
    if (KKTcheck(lambda_mid)) {
      lambda_upper <- lambda_mid
    } else {
      lambda_lower <- lambda_mid
    }
  }
  return(lambda_upper)
}
