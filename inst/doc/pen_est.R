## ----fig.show='hide'----------------------------------------------------------
library(grasps) ## for penalty computation
library(ggplot2) ## for visualization

penalties <- c("atan", "exp", "lasso", "lq", "lsp", "mcp", "scad")

pen_df <- compute_penalty(seq(-4, 4, by = 0.01), penalties, lambda = 1)
plot(pen_df, xlim = c(-1, 1), ylim = c(0, 1), zoom.size = 1) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))


## ----echo=FALSE, fig.show='hide'----------------------------------------------
plot <- plot(pen_df, xlim = c(-1, 1), ylim = c(0, 1), zoom.size = 1) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))


## ----echo=FALSE---------------------------------------------------------------
print(plot)


## ----fig.show='hide'----------------------------------------------------------
deriv_df <- compute_derivative(seq(0, 4, by = 0.01), penalties, lambda = 1)
plot(deriv_df) +
  scale_y_continuous(limits = c(0, 1.5)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))


## ----echo=FALSE, fig.show='hide'----------------------------------------------
plot <- plot(deriv_df) +
  scale_y_continuous(limits = c(0, 1.5)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))


## ----echo=FALSE---------------------------------------------------------------
print(plot)

