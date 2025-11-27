#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <algorithm> // for std::max
#include "internals.h"
using namespace Rcpp;


// Update Omega
arma::mat update_Omega(const arma::mat& M, double rho) {
  // eigenvalues decomposition
  arma::vec eigval;
  arma::mat Q;
  arma::eig_sym(eigval, Q, M);
  // update Omega
  arma::vec tildeOmega = (eigval + arma::sqrt(arma::square(eigval) + 4 * rho)) / (2 * rho);
  return Q * arma::diagmat(tildeOmega) * Q.t();
}


// Update Zblock
arma::mat update_Zblock(const arma::mat& M, const std::string& penalty,
                        const List& group_idx,
                        const double& lambda, const double& gamma,
                        const arma::mat& Initial, const double& alpha,
                        const bool diag_grp, const bool diag_include, const double rho) {

  arma::mat Z = M;
  arma::vec Z_diag = arma::diagvec(Z);
  int grp_length = group_idx.length();

  // update Z_{gh}
  for (int g = 0; g < grp_length; ++g) {
    // rows corresponding to group g
    arma::uvec row_idx = group_idx[g];
    for (int h = g; h < grp_length; ++h) {
      if (g == h && !diag_grp) {
        continue;
      }
      // columns corresponding to group h
      arma::uvec col_idx = group_idx[h];
      // get subview for in-place update
      auto Z_block = Z.submat(row_idx, col_idx);
      arma::mat Initial_block = Initial.submat(row_idx, col_idx);
      double lambda_grp = lambda_group(Initial_block, penalty, lambda, gamma);
      if (g != h) {
        double block_norm = arma::norm(Z_block, "fro");
        double scale_grp = std::max(1.0 - lambda_grp * (1-alpha) / (rho * (block_norm + eps)), 0.0);
        Z_block *= scale_grp;
        Z.submat(col_idx, row_idx) = Z_block.t();
        continue;
      }
      if (diag_include) {
        double block_norm = arma::norm(Z_block, "fro");
        double scale_grp = std::max(1.0 - lambda_grp * (1-alpha) / (rho * (block_norm + eps)), 0.0);
        Z_block *= scale_grp;
      } else {
        double sq_all = arma::accu(arma::square(Z_block));
        arma::vec de = Z_diag.elem(row_idx);
        double sq_diag = arma::dot(de, de);
        double offblock_norm = std::sqrt(std::max(sq_all-sq_diag, 0.0));
        double scale_grp = std::max(1.0 - lambda_grp * (1-alpha) / (rho * (offblock_norm + eps)), 0.0);
        Z_block *= scale_grp;
        for (arma::uword k = 0; k < row_idx.n_elem; ++k) {
          Z(row_idx[k], row_idx[k]) = de[k];
        }
      }
    }
  }

  return Z;
}

