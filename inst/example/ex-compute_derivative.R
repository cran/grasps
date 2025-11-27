library(grasps)
library(ggplot2)

deriv_df <- compute_derivative(
  omega = seq(0, 4, by = 0.01),
  penalty = c("atan", "exp", "lasso", "lq", "lsp", "mcp", "scad"),
  lambda = 1)

plot(deriv_df) +
  scale_y_continuous(limits = c(0, 1.5)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))
