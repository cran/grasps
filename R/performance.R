#' Performance Measures for Precision Matrix Estimation
#'
#' @description
#' Compute a collection of loss-based and structure-based measures to evaluate
#' the performance of an estimated precision matrix.
#'
#' @param hatOmega A numeric \eqn{d \times d} matrix giving the estimated
#' precision matrix.
#'
#' @param Omega A numeric \eqn{d \times d} matrix giving the reference
#' (typically true) precision matrix.
#'
#' @return
#' A data frame of S3 class \code{"performance"}, with one row per performance
#' metric and two columns:
#' \describe{
#' \item{measure}{The name of each performance metric. The reported metrics
#' include: sparsity, Frobenius norm loss, Kullback-Leibler divergence,
#' quadratic norm loss, spectral norm loss, true positive, true negative,
#' false positive, false negative, true positive rate, false positive rate,
#' F1 score, and Matthews correlation coefficient.}
#' \item{value}{The corresponding numeric value.}
#' }
#'
#' @details
#' Let \eqn{\Omega_{d \times d}} and \eqn{\hat{\Omega}_{d \times d}} be
#' the reference (true) and estimated precision matrices, respectively, with
#' \eqn{\Sigma = \Omega^{-1}} being the corresponding covariance matrix.
#' Edges are defined by nonzero off-diagonal entries in the upper triangle of
#' the precision matrices.
#'
#' Sparsity is treated as a structural summary, while the remaining measures
#' are grouped into loss-based measures, raw confusion-matrix counts, and
#' classification-based (structure-recovery) measures.
#'
#' "sparsity": Sparsity is computed as the proportion of zero entries among
#' the off-diagonal elements in the upper triangle of \eqn{\hat{\Omega}}.
#'
#' \strong{Loss-based measures:}
#' \itemize{
#' \item "Frobenius": Frobenius (Hilbert-Schmidt) norm loss
#' \eqn{= \Vert \Omega - \hat{\Omega} \Vert_F}.
#' \item "KL": Kullback-Leibler divergence
#' \eqn{= \mathrm{tr}(\Sigma \hat{\Omega}) - \log\det(\Sigma \hat{\Omega}) - d}.
#' \item "quadratic": Quadratic norm loss
#' \eqn{= \Vert \Sigma \hat{\Omega} - I_d \Vert_F^2}.
#' \item "spectral": Spectral (operator) norm loss
#' \eqn{= \Vert \Omega - \hat{\Omega} \Vert_{2,2} = e_1},
#' where \eqn{e_1^2} is the largest eigenvalue of \eqn{(\Omega - \hat{\Omega})^2}.
#' }
#'
#' \strong{Confusion-matrix counts:}
#' \itemize{
#' \item "TP": True positive \eqn{=} number of edges
#' in both \eqn{\Omega} and \eqn{\hat{\Omega}}.
#' \item "TN": True negative \eqn{=} number of edges
#' in neither \eqn{\Omega} nor \eqn{\hat{\Omega}}.
#' \item "FP": False positive \eqn{=} number of edges
#' in \eqn{\hat{\Omega}} but not in \eqn{\Omega}.
#' \item "FN": False negative \eqn{=} number of edges
#' in \eqn{\Omega} but not in \eqn{\hat{\Omega}}.
#' }
#'
#' \strong{Classification-based (structure-recovery) measures:}
#' \itemize{
#' \item "TPR": True positive rate (TPR), recall, sensitivity
#' \eqn{= \mathrm{TP} / (\mathrm{TP} + \mathrm{FN})}.
#' \item "FPR": False positive rate (FPR)
#' \eqn{= \mathrm{FP} / (\mathrm{FP} + \mathrm{TN})}.
#' \item "F1": \eqn{F_1} score
#' \eqn{= 2\,\mathrm{TP} / (2\,\mathrm{TP} + \mathrm{FN} + \mathrm{FP})}
#' \item "MCC": Matthews correlation coefficient (MCC)
#' \eqn{= (\mathrm{TP}\times\mathrm{TN} - \mathrm{FP}\times\mathrm{FN}) /
#' \sqrt{(\mathrm{TP}+\mathrm{FP})(\mathrm{TP}+\mathrm{FN})
#' (\mathrm{TN}+\mathrm{FP})(\mathrm{TN}+\mathrm{FN})}}
#' }
#'
#' The following table summarizes the confusion matrix and related rates:
#' \tabular{lllll}{
#' \tab \strong{Predicted Positive} \tab \strong{Predicted Negative} \tab \tab \cr
#' \strong{Real Positive} (P) \tab True positive (TP) \tab False negative (FN)
#' \tab True positive rate (TPR), recall, sensitivity = TP / P = 1 - FNR
#' \tab False negative rate (FNR) = FN / P = 1 - TPR \cr
#' \strong{Real Negative} (N) \tab False positive (FP) \tab True negative (TN)
#' \tab False positive rate (FPR) = FP / N = 1 - TNR
#' \tab True negative rate (TNR), specificity = TN / N = 1 - FPR \cr
#' \tab Positive predictive value (PPV), precision = TP / (TP + FP) = 1 - FDR
#' \tab False omission rate (FOR) = FN / (TN + FN) = 1 - NPV \tab \tab \cr
#' \tab False discovery rate (FDR) = FP / (TP + FP) = 1 - PPV
#' \tab Negative predictive value (NPV) = TN / (TN + FN) = 1 - FOR \tab \tab \cr
#' }
#'
#' @example
#' inst/example/ex-performance.R
#'
#' @export

performance <- function(hatOmega, Omega) {

  ## true covariance
  Sigma <- solve(Omega)

  ## dimension
  d <- ncol(hatOmega)

  ## Frobenius norm loss
  FL <- norm(Omega - hatOmega, "F")

  ## Kullback-Leibler divergence
  KL <- sum(diag(Sigma %*% hatOmega)) -
    determinant(Sigma %*% hatOmega, logarithm = TRUE)$modulus[1] - d

  ## quadratic loss
  QL <- norm(Sigma %*% hatOmega - diag(d), "F")^2
  ## Alternative versions:
  ## QL <- (sum(diag(Sigma %*% hatOmega - diag(d))))^2
  ## Kuismin, M. O., Kemppainen, J. T., & Sillanpää, M. J. (2017).
  ## Precision Matrix Estimation With ROPE.
  ## Journal of Computational and Graphical Statistics, 26(3), 682–694.
  ## https://doi.org/10.1080/10618600.2016.1278002
  ## QL <- sum((hatOmega %*% Sigma - diag(d))^2) ## rags2ridges

  ## spectral norm loss
  SL <- svd(Omega - hatOmega)$d[1]
  ## Equivalent version:
  ## temp <- (Omega - hatOmega) %*% (Omega - hatOmega)
  ## e1 <- max(Re(eigen(temp)$values))
  ## SL <- sqrt(e1)

  ## off-diagonal upper triangle entries define edges
  Omega_edge <- Omega[upper.tri(Omega)]
  hatOmega_edge <- hatOmega[upper.tri(hatOmega)]

  ## sparsity
  sparsity <- sum(hatOmega_edge == 0) /length(hatOmega_edge)

  ## confusion-matrix counts
  TP <- sum(Omega_edge != 0 & hatOmega_edge != 0)
  TN <- sum(Omega_edge == 0 & hatOmega_edge == 0)
  FP <- sum(Omega_edge == 0 & hatOmega_edge != 0)
  FN <- sum(Omega_edge != 0 & hatOmega_edge == 0)

  ## classification-based (structure-recovery) measures
  TPR <- TP / (TP + FN)
  ## TPR <- TP / sum(Omega_edge != 0)
  FPR <- FP / (TN + FP)
  ## FPR <- FP / sum(Omega_edge == 0)
  F1 <- 2 * TP / (2*TP + FN + FP)
  ## F1 <- 2 * precision * recall / (precision + recall)
  MCC <- (TP * TN - FP * FN) / (sqrt(TP + FP) * sqrt(TP + FN) * sqrt(TN + FP) * sqrt(TN + FN))
  ## MCC <- (TP * TN - FP * FN) / sqrt(
  ##   sum(hatOmega_edge != 0) * sum(Omega_edge != 0) * sum(Omega_edge == 0) * sum(hatOmega_edge == 0)
  ## )

  result <- data.frame(
    measure = c("sparsity", "Frobenius", "KL", "quadratic", "spectral",
                "TP", "TN", "FP", "FN", "TPR", "FPR", "F1", "MCC"),
    value = c(sparsity, FL, KL, QL, SL,
              TP, TN, FP, FN, TPR, FPR, F1, MCC)
  )
  class(result) <- c("performance", class(result))
  return(result)
}


#' @export
#' @noRd

print.performance <- function(x, ...) {
  x$value <- round(x$value, 4)
  NextMethod("print", x)
}
