#' Block-Structured Precision Matrix based on SBM
#'
#' @description
#' Generate a precision matrix that exhibits block structure induced by
#' a stochastic block model (SBM).
#'
#' @param d An integer specifying the number of variables (dimensions).
#'
#' @param block.sizes An integer vector (default = \code{NULL}) specifying
#' the size of each group. If \code{NULL}, the \eqn{d} variables are divided
#' as evenly as possible across \eqn{K} groups.
#'
#' @param K An integer (default = 3) specifying the number of groups.
#' Ignored if \code{block.sizes} is provided; then \code{K <- length(block.sizes)}.
#'
#' @param prob.mat A \eqn{K \times K} symmetric matrix (default = \code{NULL})
#' specifying the Bernoulli rates. Element \eqn{(i, j)} gives the probability of
#' creating an edge between vertices from groups \eqn{i} and \eqn{j}.
#' If \code{NULL}, a matrix with \code{within.prob} on the diagonal and
#' \code{between.prob} on the off-diagonal is used.
#'
#' @param within.prob A numeric value in [0, 1] (default = 0.25) specifying
#' the probability of creating an edge between vertices within the same group.
#' This argument is used only when \code{prob.mat = NULL}.
#'
#' @param between.prob A numeric value in [0, 1] (default = 0.05) specifying
#' the probability of creating an edge between vertices from different groups.
#' This argument is used only when \code{prob.mat = NULL}.
#'
#' @param weight.mat A \eqn{d \times d} symmetric matrix (default = \code{NULL})
#' specifying the edge weights. If \code{NULL}, weights are generated block-wise
#' according to \code{weight.dists} and \code{weight.paras}.
#'
#' @param weight.dists A list (default = \code{list("gamma", "unif")})
#' specifying the sampling distribution for each block of weights.
#' Its length determines how the distributions are assigned:
#' \itemize{
#' \item length = 1: Same specification for all blocks.
#' \item length = 2: First for within-group blocks, second for between-group
#' blocks.
#' \item length = \eqn{K + K(K-1)/2}: Full specification for each block.
#' The first \eqn{K} elements correspond to within-group blocks with indices
#' \eqn{1, \dots, K}, and the remaining \eqn{K(K-1)/2} elements correspond to
#' between-group blocks ordered as \eqn{(1,2)}, \eqn{(1,3)}, \eqn{(1,4)}, \dots,
#' \eqn{(1,K)}, \eqn{(2,3)}, \dots, \eqn{(K-1,K)}.
#' }
#' Each element of \code{weight.dists} can be:
#' \enumerate{
#' \item A user-supplied sampling function. The function must accept an argument
#' \code{n} specifying the number of samples.
#' \item A character string specifying the distribution family.
#' Accepted distributions (base R samplers in parentheses) include:
#' \itemize{
#' \item "beta": Beta distribution (\code{\link[stats]{rbeta}})
#' \item "cauchy": Cauchy distribution (\code{\link[stats]{rcauchy}}).
#' \item "chisq": Chi-squared distribution (\code{\link[stats]{rchisq}}).
#' \item "exp": Exponential distribution (\code{\link[stats]{rexp}}).
#' \item "f": F distribution (\code{\link[stats]{rf}}).
#' \item "gamma": Gamma distribution (\code{\link[stats]{rgamma}}).
#' \item "lnorm": Log normal distribution (\code{\link[stats]{rlnorm}}).
#' \item "norm": Normal distribution (\code{\link[stats]{rnorm}}).
#' \item "t": Student's t distribution (\code{\link[stats]{rt}}).
#' \item "unif": Uniform distribution (\code{\link[stats]{runif}}).
#' \item "weibull": Weibull distribution (\code{\link[stats]{rweibull}}).
#' }
#' }
#'
#' @param weight.paras A list
#' (default = \code{list(c(shape = 100, rate = 10), c(min = 0, max = 5))})
#' specifying the parameters associated with \code{weight.dists}. It must follow
#' the same length rules as \code{weight.dists}. Each element should be a named
#' vector or list suitable for the corresponding sampler.
#'
#' @param cond.target A numeric value > 1 (default = 100) specifying the target
#' condition number for the precision matrix. When necessary, a diagonal shift
#' is applied to ensure positive definiteness and numerical stability.
#'
#' @return
#' An object with S3 class \code{"gen_prec_sbm"} containing the following
#' components:
#' \describe{
#' \item{Omega}{The precision matrix with SBM block structure.}
#' \item{Sigma}{The covariance matrix, i.e., the inverse of \code{Omega}.}
#' \item{sparsity}{Proportion of zero entries in \code{Omega}.}
#' \item{membership}{An integer vector specifying the group membership.}
#' }
#'
#' @details
#' \strong{Edge sampling.}
#' Within- and between-group edges are sampled independently according to
#' Bernoulli distributions specified by \code{prob.mat}, or by \code{within.prob}
#' and \code{between.prob} if \code{prob.mat} is not supplied.
#'
#' \strong{Weight sampling.}
#' Conditional on the adjacency structure, edge weights are sampled block-wise
#' from samplers specified in \code{weight.dists} and \code{weight.paras}.
#' The length of \code{weight.dists} (and \code{weight.paras}) determines how
#' weight distributions are assigned:
#' \itemize{
#' \item length = 1: Same specification for all blocks.
#' \item length = 2: first for within-group blocks, second for between-group
#' blocks.
#' \item length = \eqn{K + K(K-1)/2}: Full specification for each block.
#' }
#'
#' \strong{Block indexing.}
#' The order for blocks is:
#' \itemize{
#' \item Within-group blocks: Indices \eqn{1, \dots, K}.
#' \item Between-group blocks: \eqn{K(K-1)/2} blocks in order \eqn{(1,2)},
#' \eqn{(1,3)}, \eqn{(1,4)}, \dots, \eqn{(1,K)}, \eqn{(2,3)}, \dots, \eqn{(K-1,K)}.
#' }
#'
#' \strong{Positive definiteness.}
#' The weighted adjacency matrix is symmetrized and used as the precision matrix
#' \eqn{\Omega_0}. Since arbitrary block-structured weights may not be positive
#' definite, a diagonal adjustment is applied to control the eigenvalue spectrum.
#' Specifically, let \eqn{\lambda_{\max}} and \eqn{\lambda_{\min}} denote
#' the largest and smallest eigenvalues of a matrix. A non-negative numeric
#' value \eqn{\tau} is added to the diagonal so that
#' \deqn{
#' \left\{
#' \begin{array}{l}
#' \dfrac{\lambda_{\max}(\Omega_0 + \tau I)}{\lambda_{\min}(\Omega_0 + \tau I)}
#' \leq \texttt{cond.target} \\[1em]
#' \lambda_{\min}(\Omega_0 + \tau I) > 0 \\[.5em]
#' \tau \geq 0
#' \end{array}
#' \right.
#' }
#' which ensures both positive definiteness and guarantees that the condition
#' number does not exceed \code{cond.target}, providing numerical stability
#' even in high-dimensional settings.
#'
#' @example
#' inst/example/ex-gen_prec_sbm.R
#'
#' @importFrom igraph as_adjacency_matrix sample_sbm
#'
#' @export

gen_prec_sbm <- function(d,
                         block.sizes = NULL, K = 3,
                         prob.mat = NULL, within.prob = 0.25, between.prob = 0.05,
                         weight.mat = NULL,
                         weight.dists = list("gamma", "unif"),
                         weight.paras = list(c(shape = 100, rate = 10),
                                             c(min = 0, max = 5)),
                         cond.target = 100) {

  ## block sizes (allow d not divisible by K)
  if (is.null(block.sizes)) {
    block.sizes <- rep(floor(d / K), K)
    rem <- d - sum(block.sizes)
    if (rem > 0) {
      block.sizes[seq_len(rem)] <- block.sizes[seq_len(rem)] + 1
    }
  }
  stopifnot(sum(block.sizes) == d)
  K <- length(block.sizes)

  ## membership and indices
  membership <- rep(1:K, times = block.sizes)
  idx_list <- split(1:d, membership)

  ## probability matrix
  if (is.null(prob.mat)) {
    prob.mat <- matrix(between.prob, K, K)
    diag(prob.mat) <- within.prob
  } else {
    stopifnot(all(dim(prob.mat) == c(K, K)))
  }

  ## SBM adjacency (undirected, no loops)
  SBM_sim <- igraph::sample_sbm(n = d, pref.matrix = prob.mat,
                                block.sizes = block.sizes,
                                directed = FALSE, loops = FALSE)
  SBM_adj <- as.matrix(igraph::as_adjacency_matrix(SBM_sim, sparse = FALSE))
  # ## generate SBM adjacency (upper triangle once, then mirror)
  # # SBM_adj <- matrix(0, d, d)
  # for (i in 1:K) {
  #   rows <- idx_list[[i]]
  #   for (j in i:K) {
  #     cols <- idx_list[[j]]
  #     block <- matrix(rbinom(length(rows) * length(cols), 1, prob.mat[i,j]),
  #                     length(rows), length(cols))
  #     SBM_adj[rows, cols] <- block
  #     if (j != i) {
  #       SBM_adj[cols, rows] <- t(block)
  #     }
  #   }
  # }
  # diag(SBM_adj) <- 0

  ## build weights
  ## fill upper-triangular; mirror to keep symmetric.
  if (is.null(weight.mat)) {
    dists_expand <- expand_spec(weight.dists, K)
    paras_expand <- expand_spec(weight.paras, K)
    weight.mat <- matrix(0, d, d)
    for (i in 1:K) {
      rows <- idx_list[[i]]
      nr <- length(rows)
      weight.mat[rows, rows] <- draw_sample(dists_expand[[i]], paras_expand[[i]], nr*nr)
      if (i < K) {
        for (j in (i+1):K) {
          cols <- idx_list[[j]]
          nc <- length(cols)
          k <- K + ((i-1)*(2*K-i))/2 + (j-i)
          weight.mat[rows, cols] <- draw_sample(dists_expand[[k]], paras_expand[[k]], nr*nc)
          weight.mat[cols, rows] <- t(weight.mat[rows, cols])
        }
      }
    }
  }

  ## precision matrix
  Omega <- SBM_adj * weight.mat
  Omega <- (Omega + t(Omega)) / 2 ## symmetric; diag still 0

  ## ensure positive definiteness and control the condition number
  ## add a scalar tau to the diagonal so that
  ## lambda_max(Omega + tau*I) / lambda_min(Omega + tau*I) <= cond.target
  ## (diagonal loading shifts all eigenvalues by the same tau)
  eigvals <- eigen(Omega, only.values = TRUE)$values
  eigval_max <- max(eigvals)
  eigval_min <- min(eigvals)
  tau <- max(0, ## no adjustment
             -eigval_min, ## ensure positive definiteness
             (eigval_max - cond.target * eigval_min) / (cond.target - 1)) ## enforce kappa <= cond.target
  diag(Omega) <- diag(Omega) + tau

  ## covariance matrix
  Sigma <- solve(Omega)

  result <- list(Omega = Omega, Sigma = Sigma,
                 sparsity = sum(Omega == 0) / length(Omega),
                 membership = membership)
  class(result) <- c("gen_prec_sbm", "blkmat")
  return(result)
}
