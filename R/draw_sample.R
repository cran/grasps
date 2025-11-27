#' Draw Random Samples from a Specified Sampler (Random Variate Generator Wrapper)
#'
#' @description
#' Generate random samples from either a base R distribution or a user-supplied
#' sampling function.
#'
#' @param sampler A specification of the sampling distribution, either:
#' \enumerate{
#' \item A user-supplied sampling function. The function must accept an argument
#' \code{n} specifying the number of samples.
#' \item A character string specifying the distribution family.
#' Supported distribution families (corresponding base R generators
#' in parentheses) include:
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
#' @param para A list (default = \code{list()}) specifying the required
#' additional parameters passed to \code{sampler}.
#'
#' @param n An integer specifying the number of samples to generate.
#'
#' @return
#' A numeric vector of length \eqn{n} containing simulated data.
#'
#' @noRd

draw_sample <- function(sampler, para = list(), n) {

  ## construct argument list
  args <- c(list(n = n), para)

  if (is.function(sampler)) { ## user-defined sampler
    tryCatch(
      do.call(sampler, args),
      error = function(e) {
        stop("Error in custom sampler: ", conditionMessage(e))
      }
    )
  } else if (is.character(sampler)) { ## base R distribution
    ## match distribution keyword to base R generator (e.g., "norm" -> "rnorm")
    # rfun <- match.fun(rfun_name)
    rfun_name <- paste0("r", sampler)
    rfun <- get0(rfun_name, mode = "function", inherits = TRUE)
    if (is.null(rfun)) {
      stop(sprintf("Unknown sampler '%s' (expected base R generator '%s').",
                   sampler, rfun_name))
    }
    tryCatch(
      do.call(rfun, args),
      error = function(e) {
        stop("Error calling ", rfun_name, ": ", conditionMessage(e))
      }
    )
  } else {
    stop("`sampler` must be either a string or a function.")
  }
}
