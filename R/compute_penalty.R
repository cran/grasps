#' Penalty Function Computation
#'
#' @description
#' Compute one or more penalty values for a given \code{omega}, allowing
#' vectorized specifications of \code{penalty}, \code{lambda}, and \code{gamma}.
#'
#' @param omega A numeric value or vector at which the penalty is evaluated.
#'
#' @param penalty A character string or vector specifying one or more penalty
#' types. Available options include:
#' \enumerate{
#' \item "lasso": Least absolute shrinkage and selection operator
#' \insertCite{tibshirani1996regression,friedman2008sparse}{grasps}.
#' \item "atan": Arctangent type penalty
#' \insertCite{wang2016variable}{grasps}.
#' \item "exp": Exponential type penalty
#' \insertCite{wang2018variable}{grasps}.
#' \item "lq": Lq penalty
#' \insertCite{frank1993statistical,fu1998penalized,fan2001variable}{grasps}.
#' \item "lsp": Log-sum penalty
#' \insertCite{candes2008enhancing}{grasps}.
#' \item "mcp": Minimax concave penalty
#' \insertCite{zhang2010nearly}{grasps}.
#' \item "scad": Smoothly clipped absolute deviation
#' \insertCite{fan2001variable,fan2009network}{grasps}.
#' }
#' If \code{penalty} has length 1, it is recycled to the common length
#' determined by \code{penalty}, \code{lambda}, and \code{gamma}.
#'
#' @param lambda A non-negative numeric value or vector specifying
#' the regularization parameter.
#' If \code{lambda} has length 1, it is recycled to the common length
#' determined by \code{penalty}, \code{lambda}, and \code{gamma}.
#'
#' @param gamma A numeric value or vector specifying the additional parameter
#' for the penalty function.
#' If \code{lambda} has length 1, it is recycled to the common length
#' determined by \code{penalty}, \code{lambda}, and \code{gamma}.
#' The penalty-specific defaults are:
#' \enumerate{
#' \item "atan": 0.005
#' \item "exp": 0.01
#' \item "lq": 0.5
#' \item "lsp": 0.1
#' \item "mcp": 3
#' \item "scad": 3.7
#' }
#' For \code{"lasso"}, \code{gamma} is ignored.
#'
#' @return
#' A data frame with S3 class \code{"penalty"} containing:
#' \describe{
#' \item{omega}{The input \code{omega} values.}
#' \item{penalty}{The penalty type for each row.}
#' \item{lambda}{The regularization parameter used.}
#' \item{gamma}{The additional penalty parameter used.}
#' \item{value}{The computed penalty value.}
#' }
#' The number of rows equals
#' \code{max(length(penalty), length(lambda), length(gamma))}.
#' Any of \code{penalty}, \code{lambda}, or \code{gamma} with length 1
#' is recycled to this common length.
#'
#' @references
#' \insertAllCited{}
#'
#' @example
#' inst/example/ex-compute_penalty.R
#'
#' @export

compute_penalty <- function(omega, penalty, lambda, gamma = NA) {

  para_len <- lengths(list(penalty, lambda, gamma))
  n <- max(para_len)
  if(para_len[1] == 1) {
    penalty <- rep(penalty, n)
  }
  if (para_len[2] == 1) {
    lambda <- rep(lambda, n)
  }
  if (para_len[3] == 1) {
    gamma <- rep(gamma, n)
  }

  res <- do.call(rbind, lapply(seq_len(n), function(k) {
    pen_internal(omega = omega, penalty = penalty[k], lambda = lambda[k], gamma = gamma[k])
  }))
  class(res) <- c("penalty", "penderiv", class(res))

  return(res)
}


#' @noRd

pen_internal <- function(omega, penalty, lambda, gamma) {

  if (!(penalty %in% c("lasso", "atan", "exp", "lq", "lsp", "mcp", "scad"))) {
    stop('Error in `penalty`!\nAvailable options: "lasso", "atan", "exp", "lq", "lsp", "mcp", "scad".')
  }

  if (lambda < 0) {
    stop('The parameter `lambda` must be non-negative!')
  }

  ## default gamma by penalty
  if (missing(gamma) || is.na(gamma)) {
    gamma <- switch(penalty,
                    "atan" = 0.005, "exp"  = 0.01, "lq"   = 0.5,
                    "lsp"  = 0.1, "mcp"  = 3, "scad" = 3.7, NA)
  }

  a <- abs(omega)

  if (penalty == "atan") {
    if (gamma <= 0) {
      warning(sprintf('For "%s", typically `gamma` > 0.', penalty), call. = FALSE)
    }
    res <- lambda * (gamma + 2/pi) * atan(a / gamma)

  } else if (penalty == "exp") {
    if (gamma <= 0) {
      warning(sprintf('For "%s", typically `gamma` > 0.', penalty), call. = FALSE)
    }
    res <- lambda * (1 - exp(-a / gamma))

  } else if (penalty == "lasso") {
    gamma <- NA
    res <- lambda * a

  } else if (penalty == "lq") {
    if (gamma <= 0 || gamma >= 1) {
      warning(sprintf('For "%s", typically 0 < `gamma` < 1.', penalty), call. = FALSE)
    }
    epsilon <- 1e-10
    res <- lambda * ((a + epsilon)^gamma)
    # res <- lambda * (pmax(a, epsilon)^gamma)
    # res <- lambda * (a^gamma)

  } else if (penalty == "lsp") {
    if (gamma <= 0) {
      warning(sprintf('For "%s", typically `gamma` > 0.', penalty), call. = FALSE)
    }
    res <- lambda * log1p(a / gamma)

  } else if (penalty == "mcp") {
    if (gamma <= 1) {
      warning(sprintf('For "%s", typically `gamma` > 1.', penalty), call. = FALSE)
    }
    res <- (lambda * a - a^2/(2*gamma)) * (a <= gamma * lambda) +
      (0.5 * gamma * lambda^2) * (a > gamma * lambda)

  } else if (penalty == "scad") {
    if (gamma <= 2) {
      warning(sprintf('For "%s", typically `gamma` > 2.', penalty), call. = FALSE)
    }
    res <- lambda * a * (a <= lambda) +
      (2 * gamma * lambda * a - a^2 - lambda^2) / (2 * (gamma-1)) * (lambda < a & a <= gamma * lambda) +
      lambda^2 * (gamma+1) / 2 * (a > gamma*lambda)
  }

  return(data.frame(omega = omega, penalty = penalty, lambda = lambda, gamma = gamma, value = res))
}
