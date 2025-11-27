#' Expand Block-Level Specifications for SBM Weight Generation
#'
#' @description
#' Expand specifications for within-group and between-group block
#' weight-generation rules into a full list of length \eqn{K + K(K-1)/2}, where:
#' - The first \eqn{K} elements correspond to within-group blocks, and
#' - The remaining \eqn{K(K-1)/2} elements correspond to between-group blocks.
#'
#' @param spec A list specifying specifications. Acceptable forms:
#' \enumerate{
#' \item length = 1: Same specification for all blocks.
#' \item length = 2: First for within-group blocks, second for between-group
#' blocks.
#' \item length = \eqn{K + K(K-1)/2}: Full specification for each block.
#' }
#'
#' @param K An integer specifying the number of groups.
#'
#' @return
#' A list of length \eqn{K + K(K-1)/2}, where the first \eqn{K} elements
#' correspond to within-group blocks, and the remaining \eqn{K(K-1)/2}
#' correspond to between-group blocks.
#'
#' @noRd

expand_spec <- function(spec, K) {

  ## number of within-group and between-group blocks
  n_within <- K
  n_between <- K*(K-1)/2
  n_total <- n_within + n_between

  ## ensure the input is a list
  ## if it is a vector, convert to a list
  if (!is.list(spec)) {
    spec <- as.list(spec)
  }

  spec_len <- length(spec)
  if (spec_len == 1) {
    ## Case 1: one rule -> use the same specification for all blocks
    return(rep(list(spec[[1]]), n_total))

  } else if (spec_len == 2) {
    ## Case 2: two rules -> 1st rule for all within-group blocks
    ##                      2nd rule for all between-group blocks
    return(c(rep(list(spec[[1]]), n_within),
             rep(list(spec[[2]]), n_between)))

  } else if (spec_len == n_total) {
    ## Case 3: full specification explicitly provided
    return(spec)
  }

  ## otherwise: invalid specification length
  stop(paste0("`weight.dists`/`weight.paras` length must be 1, 2, or K + K*(K-1)/2.\n",
              sprintf("For K = %d, expected lengths: 1, 2, or %d.", K, n_total)))
}
