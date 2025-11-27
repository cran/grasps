#' Groupwise Regularized Adaptive Sparse Precision Solution
#'
#' @description
#' Provide a collection of statistical methods that incorporate both
#' element-wise and group-wise penalties to estimate a precision matrix.
#'
#' @param X \enumerate{
#' \item An \eqn{n \times d} data matrix with sample size \eqn{n} and
#' dimension \eqn{d}.
#' \item A \eqn{d \times d} sample covariance matrix with dimension \eqn{d}.
#' }
#'
#' @param n An integer (default = \code{nrow(X)}) specifying the sample size.
#' This is only required when the input matrix \code{X} is a \eqn{d \times d}
#' sample covariance matrix with dimension \eqn{d}.
#'
#' @param membership An integer vector specifying the group membership.
#' The length of \code{membership} must be consistent with the dimension \eqn{d}.
#'
#' @param penalty A character string specifying the penalty for estimating
#' precision matrix. Available options include:
#' \enumerate{
#' \item "lasso": Least absolute shrinkage and selection operator
#' \insertCite{tibshirani1996regression,friedman2008sparse}{grasps}.
#' \item "adapt": Adaptive lasso
#' \insertCite{zou2006adaptive,fan2009network}{grasps}.
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
#'
#' @param diag.ind A logical value (default = TRUE) specifying whether to
#' penalize the diagonal elements.
#'
#' @param diag.grp A logical value (default = TRUE) specifying whether to
#' penalize the within-group blocks.
#'
#' @param diag.include A logical value (default = FALSE) specifying whether to
#' include the diagonal entries in the penalty for within-group blocks when
#' \code{diag.grp = TRUE}.
#'
#' @param lambda A non-negative numeric vector specifying the grid for
#' the regularization parameter. The default is \code{NULL}, which generates
#' its own \code{lambda} sequence based on \code{nlambda} and
#' \code{lambda.min.ratio}.
#'
#' @param alpha A numeric vector in [0, 1] specifying the grid for
#' the mixing parameter balancing the element-wise individual L1 penalty and
#' the block-wise group L2 penalty.
#' An alpha of 1 corresponds to the individual penalty only; an alpha of 0
#' corresponds to the group penalty only.
#' The default value is a sequence from 0.1 to 0.9 with increments of 0.1.
#'
#' @param gamma A numeric value specifying the additional parameter fo
#' the chosen \code{penalty}. The default value depends on the penalty:
#' \enumerate{
#' \item "adapt": 0.5
#' \item "atan": 0.005
#' \item "exp": 0.01
#' \item "lq": 0.5
#' \item "lsp": 0.1
#' \item "mcp": 3
#' \item "scad": 3.7
#' }
#'
#' @param nlambda An integer (default = 10) specifying the number of
#' \code{lambda} values to generate when \code{lambda = NULL}.
#'
#' @param lambda.min.ratio A numeric value > 0 (default = 0.01) specifying
#' the fraction of the maximum \code{lambda} value \eqn{\lambda_{max}} to
#' generate the minimum \code{lambda} \eqn{\lambda_{min}}.
#' If \code{lambda = NULL}, a \code{lambda} grid of length \code{nlambda} is
#' automatically generated on a log scale, ranging from \eqn{\lambda_{max}}
#' down to \eqn{\lambda_{min}}.
#'
#' @param growiter.lambda An integer (default = 30) specifying the maximum
#' number of exponential growth steps during the initial search for an
#' admissible upper bound \eqn{\lambda_{\max}}.
#'
#' @param tol.lambda A numeric value > 0 (default = 1e-03) specifying
#' the relative tolerance for the bisection stopping rule on the interval width.
#'
#' @param maxiter.lambda An integer (default = 50) specifying the maximum number
#' of bisection iterations in the line search for  \eqn{\lambda_{\max}}.
#'
#' @param rho A numeric value > 0 (default = 2) specifying the ADMM
#' augmented-Lagrangian penalty parameter (often called the ADMM step size).
#' Larger values typically put more weight on enforcing the consensus
#' constraints at each iteration; smaller values yield more conservative updates.
#'
#' @param tau.incr A numeric value > 1 (default = 2) specifying
#' the multiplicative factor used to increase \code{rho} when the primal
#' residual dominates the dual residual in ADMM.
#'
#' @param tau.decr A numeric value > 1 (default = 2) specifying
#' the multiplicative factor used to decrease \code{rho} when the dual residual
#' dominates the primal residual in ADMM.
#'
#' @param nu A numeric value > 1 (default = 10) controlling how aggressively
#' \code{rho} is rescaled in the adaptive-\code{rho} scheme (residual balancing).
#'
#' @param tol.abs A numeric value > 0 (default = 1e-04) specifying the absolute
#' tolerance for ADMM stopping (applied to primal/dual residual norms).
#'
#' @param tol.rel A numeric value > 0 (default = 1e-04) specifying the relative
#' tolerance for ADMM stopping (applied to primal/dual residual norms).
#'
#' @param maxiter An integer (default = 1e+04) specifying the maximum number of
#' ADMM iterations.
#'
#' @param crit A character string (default = "BIC") specifying the parameter
#' selection criterion to use. Available options include:
#' \enumerate{
#' \item "AIC": Akaike information criterion
#' \insertCite{akaike1973information}{grasps}.
#' \item "BIC": Bayesian information criterion
#' \insertCite{schwarz1978estimating}{grasps}.
#' \item "EBIC": extended Bayesian information criterion
#' \insertCite{chen2008extended,foygel2010extended}{grasps}.
#' \item "HBIC": high dimensional Bayesian information criterion
#' \insertCite{wang2013calibrating,fan2017high}{grasps}.
#' \item "CV": k-fold cross validation with negative log-likelihood loss.
#' }
#'
#' @param kfold An integer (default = 5) specifying the number of folds used for
#' \code{crit = "CV"}.
#'
#' @param ebic.tuning A numeric value in [0, 1] (default = 0.5) specifying
#' the tuning parameter to calculate for \code{crit = "EBIC"}.
#'
#' @return
#' An object with S3 class \code{"grasps"} containing the following components:
#' \describe{
#' \item{hatOmega}{The estimated precision matrix.}
#' \item{lambda}{The optimal regularization parameter.}
#' \item{alpha}{The optimal mixing parameter.}
#' \item{initial}{The initial estimate of \code{hatOmega} when a non-convex
#' penalty is chosen via \code{penalty}.}
#' \item{gamma}{The optimal addtional parameter when a non-convex penalty
#' is chosen via \code{penalty}.}
#' \item{iterations}{The number of ADMM iterations.}
#' \item{lambda.grid}{The actual lambda grid used in the program.}
#' \item{alpha.grid}{The actual alpha grid used in the program.}
#' \item{lambda.safe}{The bisection-refined upper bound \eqn{\lambda_{\max}},
#' corresponding to \code{alpha.grid}, when \code{lambda = NULL}.}
#' \item{loss}{The optimal k-fold loss when \code{crit = "CV"}.}
#' \item{CV.loss}{Matrix of CV losses, with rows for parameter combinations and
#' columns for CV folds, when \code{crit = "CV"}.}
#' \item{score}{The optimal information criterion score when \code{crit} is set
#' to \code{"AIC"}, \code{"BIC"}, \code{"EBIC"}, or \code{"HBIC"}.}
#' \item{IC.score}{The information criterion score for each parameter
#' combination when \code{crit} is set to \code{"AIC"}, \code{"BIC"},
#' \code{"EBIC"}, or \code{"HBIC"}.}
#' \item{membership}{The group membership.}
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @example
#' inst/example/ex-grasps.R
#'
#' @importFrom stats cov
#' @importFrom Rdpack reprompt
#'
#' @export

grasps <- function(X, n = nrow(X), membership, penalty,
                   diag.ind = TRUE, diag.grp = TRUE, diag.include = FALSE,
                   lambda = NULL, alpha = NULL, gamma = NULL,
                   nlambda = 10, lambda.min.ratio = 0.01,
                   growiter.lambda = 30, tol.lambda = 1e-03, maxiter.lambda = 50,
                   rho = 2, tau.incr = 2, tau.decr = 2, nu = 10,
                   tol.abs = 1e-04, tol.rel = 1e-04, maxiter = 1e+04,
                   crit = "BIC", kfold = 5, ebic.tuning = 0.5) {

  d <- ncol(X)

  if (length(membership) != d) {
    stop('The length of `membership` must equal the column dimension of `X`!')
  }

  if (!(penalty %in% c("lasso", "adapt", "atan", "exp", "lq", "lsp", "mcp", "scad"))) {
    stop('Error in `penalty`!
         Available options: "lasso", "adapt", "atan", "exp", "lq", "lsp", "mcp", "scad".')
  }

  if (!(crit %in% c("AIC", "BIC", "EBIC", "HBIC", "CV"))) {
    stop('Error in `crit`!
         Available options: "AIC", "BIC", "EBIC", "HBIC", "CV".')
  }

  if (!all(lambda > 0)) {
    stop('The parameter `lambda` must be positive!')
  }

  if (!all(alpha >= 0 & alpha <= 1)) {
    stop('The parameter `alpha` must be in [0,1]!')
  }

  if (rho <= 0) {
    stop('The parameter `rho` must be positive!')
  }

  if (tau.incr <= 1) {
    stop('The parameter `tau.incr` must be greater than 1!')
  }

  if (tau.decr <= 1) {
    stop('The parameter `tau.decr` must be greater than 1!')
  }

  if (nu <= 1) {
    stop('The parameter `nu` must be greater than 1!')
  }

  if (tol.abs <= 0) {
    stop('The parameter `tol.abs` must be positive!')
  }

  if (tol.rel <= 0) {
    stop('The parameter `tol.rel` must be positive!')
  }

  if (isSymmetric(X)) {
    S <- X
    X <- NULL
  } else {
    S <- (n-1)/n*cov(X)
  }

  grp <- sort(unique(membership))
  group.idx <- lapply(grp, function(g) which(membership == g)-1)
  names(group.idx) <- grp

  alpha_null <- is.null(alpha)
  lambda_null <- is.null(lambda)

  if (alpha_null) {
    alpha <- seq(0.1, 0.9, 0.1)
  }

  if (lambda_null) {
    max_abs.off <- max(abs(S[upper.tri(S, diag = FALSE)]))
    S2 <- S*S
    sum.row <- rowsum(S2, membership)
    block.sum <- rowsum(t(sum.row), membership)
    block.norm <- sqrt(block.sum)
    max_block.norm.off <- max(block.norm[upper.tri(block.norm, diag = FALSE)])
    idx_list <- split(seq_along(membership), membership)
    block.off_list <- list()
    for (g in seq_along(idx_list)) {
      for (gp in seq_along(idx_list)) {
        if (g < gp) {
          block.off_list[[length(block.off_list)+1]] <- S[idx_list[[g]], idx_list[[gp]], drop = FALSE]
        }
      }
    }
    lambda_path <- lapply(alpha, function(a) {
      lambda.ind <- ifelse(a > 0, max_abs.off / a, 0)
      lambda.grp <- ifelse(a < 1, max_block.norm.off / (1-a), 0)
      lambda.safe <- max(lambda.ind, lambda.grp)
      if (a > 0 & a < 1) {
        lambda.max <- line_search_lambda_max(block.off_list, lambda.safe, a, growiter.lambda, tol.lambda, maxiter.lambda)
      } else {
        lambda.max <- lambda.safe
      }
      lambda.min <- lambda.min.ratio*lambda.max
      return(list(lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda)),
                  lambda.safe = lambda.safe))
    })
    lambda <- unlist(lapply(lambda_path, `[[`, "lambda"))
    lambda.safe <- sapply(lambda_path, `[[`, "lambda.safe")
    names(lambda.safe) <- alpha
  }

  if (lambda_null & alpha_null) {
    parameter <- data.frame(alpha = rep(alpha, each = nlambda),
                            lambda = lambda)
  } else {
    parameter <- data.frame(alpha = alpha, lambda = lambda)
  }

  ## default gamma by penalty
  if (missing(gamma) || is.null(gamma)) {
    gamma <- switch(penalty,
                    "adapt" = 0.5, "atan" = 0.005, "exp"  = 0.01, "lq"   = 0.5,
                    "lsp"  = 0.1, "mcp"  = 3, "scad" = 3.7, NA)
  }

  if (nrow(parameter) > 1) {

    if (crit == "CV") {

      if(is.null(X)) {
        stop('CV requires the n-by-d data matrix!')
      }

      if (kfold < 2 | kfold > n) {
        stop('`kfold` must be between 2 and the row dimension of `X`!')
      }

      CV <- ADMMsggm_CV(X = X, group_idx = group.idx, penalty = penalty,
                        diag_ind = diag.ind, diag_grp = diag.grp, diag_include = diag.include,
                        lambdas = parameter$lambda, alphas = parameter$alpha, gamma = gamma,
                        rho = rho, tau_incr = tau.incr, tau_decr = tau.decr, nu = nu,
                        tol_abs = tol.abs, tol_rel = tol.rel, maxiter = maxiter,
                        kfold = kfold)

      result <- ADMMsggm(S = S, group_idx = group.idx, penalty = penalty,
                         diag_ind = diag.ind, diag_grp = diag.grp, diag_include = diag.include,
                         lambda = CV$lambda.opt, alpha = CV$alpha.opt, gamma = gamma,
                         rho = rho, tau_incr = tau.incr, tau_decr = tau.decr, nu = nu,
                         tol_abs = tol.abs, tol_rel = tol.rel, maxiter = maxiter)
      result$lambda.grid <- parameter$lambda
      result$alpha.grid <- parameter$alpha
      if (lambda_null) {
        result$lambda.safe <- lambda.safe
      }
      result$loss <- CV$loss.min
      result$CV.loss <- CV$CV.loss

    } else {

      IC <- ADMMsggm_IC(S = S, group_idx = group.idx, penalty = penalty,
                        diag_ind = diag.ind, diag_grp = diag.grp, diag_include = diag.include,
                        lambdas = parameter$lambda, alphas = parameter$alpha, gamma = gamma,
                        rho = rho, tau_incr = tau.incr, tau_decr = tau.decr, nu = nu,
                        tol_abs = tol.abs, tol_rel = tol.rel, maxiter = maxiter,
                        crit = crit, n = n, ebic_tuning = ebic.tuning)

      result <- IC$result
      result$lambda.grid <- parameter$lambda
      result$alpha.grid <- parameter$alpha
      if (lambda_null) {
        result$lambda.safe <- lambda.safe
      }
      result$score <- IC$score.min
      result$IC.score <- as.vector(IC$IC.score)

    }

  } else {

    result <- ADMMsggm(S = S, group_idx = group.idx, penalty = penalty,
                       diag_ind = diag.ind, diag_grp = diag.grp, diag_include = diag.include,
                       lambda = lambda, alpha = alpha, gamma = gamma,
                       rho = rho, tau_incr = tau.incr, tau_decr = tau.decr, nu = nu,
                       tol_abs = tol.abs, tol_rel = tol.rel, maxiter = maxiter)

    if (lambda_null) {
      result$lambda.safe <- lambda.safe
    }

  }

  result$membership <- membership
  class(result) <- c("grasps", "blkmat")
  return(result)

}
