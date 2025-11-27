#' Plot Function for S3 Class "penderiv"
#'
#' @description
#' Generate a visualization of penalty functions produced by
#' \code{\link[grasps]{compute_penalty}}, or penalty derivatives produced by
#' \code{\link[grasps]{compute_derivative}}.
#' The plot automatically summarizes multiple configurations of penalty type,
#' \eqn{\lambda}, and \eqn{\gamma}. Optional zooming is supported through
#' \code{\link[ggforce]{facet_zoom}}.
#'
#' @param x An object inheriting from S3 class \code{"penderiv"}, typically
#' returned by \code{\link[grasps]{compute_penalty}}, or
#' \code{\link[grasps]{compute_derivative}}.
#'
#' @param ... Optional arguments passed to \code{\link[ggforce]{facet_zoom}}
#' to zoom in on a subset of the data, while keeping the view of the full
#' dataset as a separate panel.
#'
#' @return
#' An object of class \code{ggplot}.
#'
#' @example
#' inst/example/ex-plot.penderiv.R
#'
#' @import ggplot2
#' @import ggforce
#'
#' @export

plot.penderiv <- function(x, ...) {

  ## identify which configuration columns vary across rows
  target <- c("penalty", "lambda", "gamma")
  cols_keep <- target[sapply(x[target], function(col) {
    col_na_rm <- col[!is.na(col)]
    length(col_na_rm) && (anyNA(col) || any(col_na_rm != col_na_rm[1]))
  })]

  ## clean data
  df <- x[, (names(x) %in% c("omega", "value", cols_keep)), drop = FALSE]
  key <- paste(cols_keep, collapse = "_")

  ## build facet_zoom() layer if zooming arguments are provided
  fz <- if (length(list(...)) > 0L) {
    do.call(facet_zoom, list(...))
  } else {
    NULL
  }

  ## declare
  group <- omega <- value <- NULL

  ## y-axis label
  y_lab <- if (inherits(x, "penalty")) {
    expression("Penalty Function" ~ italic(lambda) * italic(p) ~ "(" * italic(omega) * ")")
  } else if (inherits(x, "derivative")) {
    expression("Derivative Function" ~ italic(lambda) * italic(p) ~ "'(" * italic(omega) * ")")
  } else {
    NULL
  }

  ## all configurations identical (only one curve)
  if (key == "") {
    ggplot(df, aes(x = omega, y = value)) +
      geom_line() +
      fz +
      labs(x = expression(italic(omega)), y = y_lab) +
      theme_bw() +
      theme(legend.position = "bottom")

  ## multiple configurations, use color grouping and legend
  } else {

    ## group labels based on which parameter(s) vary
    df$group <- switch(
      key,
      "penalty" = df$penalty,
      "lambda" = sprintf("paste(lambda, ' = ', %s)", df$lambda),
      "gamma" = sprintf("paste(gamma,  ' = ', %s)", df$gamma),
      "penalty_lambda" = sprintf(
        "paste('%s', ' (', lambda, ' = ', %s, ')')", df$penalty, df$lambda
      ),
      "penalty_gamma" = ifelse(
        df$penalty == "lasso",
        df$penalty,
        sprintf("paste('%s', ' (', gamma, ' = ', %s, ')')", df$penalty, df$gamma)
      ),
      "lambda_gamma" = ifelse(
        df$penalty == "lasso",
        sprintf("paste(lambda, ' = ', %s)", df$lambda),
        sprintf("paste(lambda, ' = ', %s, ', ', gamma, ' = ', %s)", df$lambda, df$gamma)
      ),
      "penalty_lambda_gamma" = ifelse(
        df$penalty == "lasso",
        sprintf("paste('%s', ' (', lambda, ' = ', %s, ')')", df$penalty, df$lambda),
        sprintf("paste('%s', ' (', lambda, ' = ', %s, ', ', gamma, ' = ', %s, ')')",
                df$penalty, df$lambda, df$gamma)
      )
    )

    ## legend label based on which parameter(s) vary
    legend_name <- switch(
      key,
      lambda = expression(lambda),
      gamma  = expression(gamma),
      "Penalty Type"
    )

    ggplot(df, aes(x = omega, y = value, color = group)) +
      geom_line() +
      fz +
      labs(x = expression(italic(omega)), y = y_lab, color = legend_name) +
      scale_color_discrete(labels = function(z) parse(text = z)) +
      theme_bw() +
      theme(legend.position = "bottom")
  }
}
