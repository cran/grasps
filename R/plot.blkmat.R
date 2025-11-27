#' Plot Function for Block-Structured Precision Matrices
#' (Visualize a Matrix with Group Boundaries)
#'
#' @description
#' Visualize a precision matrix as a heatmap with dashed boundary lines
#' separating group blocks. This function is shared by objects returned from
#' \code{\link[grasps]{grasps}}, \code{\link[grasps]{gen_prec_sbm}}, and
#' \code{\link[grasps]{sparsify_block_banded}}, all of which inherit from
#' the S3 class \code{"blkmat"}.
#'
#' @param x An object inheriting from S3 class \code{"blkmat"}, typically
#' returned by \code{\link[grasps]{grasps}}, \code{\link[grasps]{gen_prec_sbm}}
#' or \code{\link[grasps]{sparsify_block_banded}}.
#'
#' @param colors A vector of colors specifying an n-color gradient scale for
#' the fill aesthetics.
#'
#' @param ... Additional arguments passed to \code{\link[ggplot2]{ggplot}}.
#'
#' @return
#' A heatmap of class \code{ggplot} showing the matrix entries.
#' Dashed lines indicate group boundaries.
#' The plot title also reports matrix dimension and sparsity.
#'
#' @example
#' inst/example/ex-plot.blkmat.R
#'
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom scales rescale
#'
#' @export

plot.blkmat <- function(x, colors = NULL, ...) {

  if (is.null(colors)) {
    colors <- colorRampPalette(
      c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
        "yellow", "#FF7F00", "red", "#7F0000"))(256)
  }

  ## choose which matrix to plot
  if (inherits(x, "grasps")) {
    mat <- x$hatOmega
  } else if (inherits(x, c("gen_prec_sbm", "sparsify_block_banded"))) {
    mat <- x$Omega
  }

  d <- ncol(mat)

  ## compute group sizes and boundary positions for dashed lines
  grp_sizes <- table(x$membership)
  cuts <- cumsum(grp_sizes)
  bnds <- cuts[-length(cuts)] + 0.5 ## boundaries between groups
  y_bnds <- d - bnds + 1  ## flipped y, note: scale_y_discrete(limits = rev)

  ## declare
  Col <- Row <- value <- NULL

  ## plot data
  labs <- paste0("V", seq_len(d))
  plotData <- data.frame(
    Row = factor(rep(labs, times = d),  levels = labs),
    Col = factor(rep(labs, each = d), levels = labs),
    value = as.vector(mat),
    check.names = FALSE
  )
  # plotData <- as.data.frame(mat) %>%
  #   mutate(Row = factor(paste0("V", rownames(.)),
  #                       levels = paste0("V", rownames(.)))) %>%
  #   pivot_longer(-Row, names_to = "Col", values_to = "value") %>%
  #   mutate(Col = factor(Col, levels = levels(Row)))

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
    scale_fill_gradientn(colours = colors,
                         values = rescale(c(vmin, 0, vmax)),
                         limits = c(vmin, vmax),
                         na.value = "white") +
    geom_vline(xintercept = bnds, linetype = "dashed") +
    geom_hline(yintercept = y_bnds, linetype = "dashed") +
    labs(x = NULL, y = NULL,
         title = sprintf("Dimension = %d, Sparsity = %s", d, round(sparsity, 4))) +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.ticks = element_blank(),
          # plot.margin = margin(1, 1, 1, 1, "mm"),
          plot.title = element_text(hjust = .5))
}
