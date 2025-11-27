#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "internals.h"
using namespace Rcpp;


// Element-Wise Tuning Parameter Matrix of Individual Penalty
arma::mat lambda_individual(const arma::mat& Omega, std::string penalty,
                            double lambda, double gamma) {

  arma::mat W = arma::zeros(Omega.n_rows, Omega.n_cols);

  arma::mat abs_Omega = arma::abs(Omega);

  if (penalty == "lasso") {
    W.fill(lambda);

  } else if (penalty == "adapt") {
    W = lambda / (arma::pow(abs_Omega, gamma) + eps);

  } else if (penalty == "atan") {
    W = lambda * gamma * (gamma + 2/M_PI) / (gamma*gamma + Omega%Omega);

  } else if (penalty == "exp") {
    W = (lambda / gamma) * arma::exp(-abs_Omega / gamma);

  } else if (penalty == "lq") {
    W = lambda * gamma * arma::pow(abs_Omega + eps, gamma - 1.0);

  } else if (penalty == "lsp") {
    W = lambda / (gamma + abs_Omega);

  } else if (penalty == "mcp") {
    W = (lambda - abs_Omega/gamma) % (abs_Omega <= gamma*lambda);

  } else if (penalty == "scad") {
    W.fill(lambda);
    W = W % (abs_Omega <= lambda);
    W += ((gamma*lambda - abs_Omega) / (gamma - 1)) %
      ((abs_Omega > lambda) % (abs_Omega <= gamma*lambda));

  } else {
    stop("Unsupported penalty type!");
  }

  return W;
}


// Tuning Parameter Scalar of Group Penalty
double lambda_group(const arma::mat& Omega, std::string penalty,
                    double lambda, double gamma) {

  double Omega_norm = arma::norm(Omega, "fro");

  if (penalty == "lasso") {
    return lambda;

  } else if (penalty == "adapt") {
    return lambda / (arma::norm(arma::pow(arma::abs(Omega), gamma), "fro") + eps);

  } else if (penalty == "atan") {
    return lambda * gamma * (gamma + 2/M_PI) / (gamma*gamma + Omega_norm*Omega_norm);

  } else if (penalty == "exp") {
    return (lambda / gamma) * std::exp(-Omega_norm / gamma);

  } else if (penalty == "lq") {
    return lambda * gamma * std::pow(Omega_norm + eps, gamma - 1.0);

  } else if (penalty == "lsp") {
    return lambda / (gamma + Omega_norm);

  } else if (penalty == "mcp") {
    if (Omega_norm <= gamma*lambda) {
      return lambda - Omega_norm/gamma;
    } else {
      return 0.0;
    }

  } else if (penalty == "scad") {
    if (Omega_norm <= lambda) {
      return lambda;
    } else if (Omega_norm <= gamma*lambda) {
      return (gamma*lambda - Omega_norm) / (gamma - 1);
    } else {
      return 0.0;
    }

  } else {
    stop("Unsupported penalty type!");
  }
}

