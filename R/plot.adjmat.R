#' Plot Function for S3 Class "adjmat"
#'
#' @description
#' Visualize an adjacency matrix as a heatmap. This function is shared by
#' objects returned from \code{\link[grasps]{prec_to_adj}}.
#'
#' @param x An object inheriting from S3 class \code{"adjmat"}, typically
#' returned by \code{\link[grasps]{prec_to_adj}}.
#'
#' @param ... Additional arguments passed to \code{\link[ggplot2]{ggplot}}.
#'
#' @return
#' A heatmap of class \code{ggplot} showing the matrix entries.
#' The plot title also reports matrix dimension and sparsity.
#'
#' @example
#' inst/example/ex-grasps.R
#'
#' @import ggplot2
#'
#' @export

plot.adjmat <- function(x, ...) {

  mat <- x
  p <- ncol(mat)

  ## declare
  Col <- Row <- value <- NULL

  ## plot data
  labs <- as.character(seq_len(p))
  plotData <- data.frame(
    Row = factor(rep(labs, times = p),  levels = labs),
    Col = factor(rep(labs, each = p), levels = labs),
    value = as.vector(mat),
    check.names = FALSE
  )

  ## zero -> NA for better white rendering
  plotData$value[plotData$value == 0] <- NA

  ## color scaling range
  vmin <- min(plotData$value, na.rm = TRUE)
  vmax <- max(plotData$value, na.rm = TRUE)

  ## sparsity
  mat_edge <- mat[upper.tri(mat)]
  sparsity <- sum(mat_edge == 0) /length(mat_edge)

  ggplot(plotData, aes(x = Col, y = Row, fill = value)) +
    coord_fixed() +
    geom_tile() +
    guides(fill = guide_colourbar(title = NULL, barwidth = 0.5, barheight = 5)) +
    scale_x_discrete(limits = labs, expand = c(0, 0)) +
    scale_y_discrete(limits = rev(labs), expand = c(0, 0)) +
    scale_fill_viridis_c(option = "turbo",
                         values = rescale(c(vmin, 0, vmax)),
                         limits = c(vmin, vmax),
                         na.value = "white") +
    labs(x = NULL, y = NULL,
         title = sprintf("Dimension = %d, Sparsity = %.4f", p, sparsity)) +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.ticks = element_blank(),
          # plot.margin = margin(1, 1, 1, 1, "mm"),
          plot.title = element_text(hjust = .5))
}
