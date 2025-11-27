#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <algorithm> // for std::max
#include <cmath> // for std::isfinite()
#include "internals.h"
using namespace Rcpp;


// Sparse-Group Graphical Lasso via ADMM
List ADMMsggl(const arma::mat& S, const List& group_idx,
              bool diag_ind, bool diag_grp, bool diag_include,
              const double& lambda, const double& alpha,
              double rho, const double tau_incr, const double tau_decr, const double nu,
              const double tol_abs, const double tol_rel, const int maxiter) {

  int d = S.n_cols;
  int iter = 0;
  arma::mat Omega = arma::diagmat(1/arma::diagvec(S));
  arma::mat Z = arma::zeros(d, d);
  arma::mat U = arma::zeros(d, d);
  arma::mat Z_old = Z;

  arma::mat I = arma::ones(d, d);
  if (!diag_ind) {
    I -= arma::diagmat(I);
  }

  int grp_length = group_idx.length();

  // loop until convergence
  for (iter = 0; iter < maxiter; ++iter) {

    Z_old = Z;

    // update Omega
    Omega = update_Omega(rho*(Z-U)-S, rho);

    // update Z
    // step 1: element-wise individual L1 penalty
    Z = soft_matrix(Omega + U, lambda * I * alpha / rho);
    arma::vec Z_diag = arma::diagvec(Z);
    // step 2: block-wise group L2 penalty
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
        if (g != h) {
          double block_norm = arma::norm(Z_block, "fro");
          double scale_grp = std::max(1.0 - lambda * (1-alpha) / (rho * (block_norm + eps)), 0.0);
          Z_block *= scale_grp;
          Z.submat(col_idx, row_idx) = Z_block.t();
          continue;
        }
        if (diag_include) {
          double block_norm = arma::norm(Z_block, "fro");
          double scale_grp = std::max(1.0 - lambda * (1-alpha) / (rho * (block_norm + eps)), 0.0);
          Z_block *= scale_grp;
        } else {
          double sq_all = arma::accu(arma::square(Z_block));
          arma::vec de = Z_diag.elem(row_idx);
          double sq_diag = arma::dot(de, de);
          double offblock_norm = std::sqrt(std::max(sq_all-sq_diag, 0.0));
          double scale_grp = std::max(1.0 - lambda * (1-alpha) / (rho * (offblock_norm + eps)), 0.0);
          Z_block *= scale_grp;
          for (arma::uword k = 0; k < row_idx.n_elem; ++k) {
            Z(row_idx[k], row_idx[k]) = de[k];
          }
        }
      }
    }

    // update U
    U += Omega - Z;

    // residual
    double r = arma::norm(Omega - Z, "fro");
    double s = arma::norm(rho*(Z - Z_old), "fro");

    // tolerance
    double tol_pri = d * tol_abs + tol_rel * std::max(arma::norm(Omega, "fro"), arma::norm(Z, "fro"));
    double tol_dual = d * tol_abs + tol_rel * rho * arma::norm(U, "fro");

    if (r <= tol_pri && s <= tol_dual) {
      break;
    }

    // update rho
    if (r > nu * s) {
      rho *= tau_incr;
      U /= tau_incr;
    }
    if (s > nu * r) {
      rho /= tau_decr;
      U *= tau_decr;
    }

    // R_CheckUserInterrupt
    if (iter % 1000 == 0) {
      R_CheckUserInterrupt();
    }
  }

  return List::create(
    Named("hatOmega") = Z,
    Named("lambda") = lambda,
    Named("alpha") = alpha,
    Named("iterations") = iter
  );
}


// Sparse-Group Graphical Non-Convex Penalties via ADMM
List ADMMsggn(const arma::mat& S, const List& group_idx, std::string penalty,
              bool diag_ind, bool diag_grp, bool diag_include,
              const double& lambda, const double& alpha, const double& gamma,
              double rho, const double tau_incr, const double tau_decr, const double nu,
              const double tol_abs, const double tol_rel, const int maxiter) {

  int d = S.n_cols;
  int iter = 0;
  arma::mat Omega = arma::diagmat(1/arma::diagvec(S));
  arma::mat Z = arma::zeros(d, d);
  arma::mat U = arma::zeros(d, d);
  arma::mat Z_old = Z;
  List sgglres = ADMMsggl(S, group_idx, diag_ind, diag_grp, diag_include, lambda, alpha,
                          rho, tau_incr, tau_decr, nu,
                          tol_abs, tol_rel, maxiter);
  arma::mat initial = sgglres["hatOmega"];

  arma::mat I = arma::ones(d, d);
  if (!diag_ind) {
    I -= arma::diagmat(I);
  }
  arma::mat lambda_ind = lambda_individual(initial, penalty, lambda, gamma) % I;

  // loop until convergence
  for (iter = 0; iter < maxiter; ++iter) {

    Z_old = Z;

    // update Omega
    Omega = update_Omega(rho*(Z-U)-S, rho);

    // update Z
    // step 1: element-wise individual L1 penalty
    Z = soft_matrix(Omega + U, lambda_ind * alpha / rho);

    // step 2: block-wise group L2 penalty
    Z = update_Zblock(Z, penalty, group_idx, lambda, gamma, initial, alpha, diag_grp, diag_include, rho);

    // update U
    U += Omega - Z;

    // residual
    double r = arma::norm(Omega - Z, "fro");
    double s = arma::norm(rho*(Z - Z_old), "fro");

    // tolerance
    double tol_pri = d * tol_abs + tol_rel * std::max(arma::norm(Omega, "fro"), arma::norm(Z, "fro"));
    double tol_dual = d * tol_abs + tol_rel * rho * arma::norm(U, "fro");

    if (r <= tol_pri && s <= tol_dual) {
      break;
    }

    // update rho
    if (r > nu * s) {
      rho *= tau_incr;
      U /= tau_incr;
    }
    if (s > nu * r) {
      rho /= tau_decr;
      U *= tau_decr;
    }

    // R_CheckUserInterrupt
    if (iter % 1000 == 0) {
      R_CheckUserInterrupt();
    }
  }

  return List::create(
    Named("hatOmega") = Z,
    Named("lambda") = lambda,
    Named("alpha") = alpha,
    Named("initial") = initial,
    Named("gamma") = gamma,
    Named("iterations") = iter
  );
}


// Sparse-Group Graphical Model via ADMM
// [[Rcpp::export]]
List ADMMsggm(const arma::mat& S, const List& group_idx, std::string penalty,
              bool diag_ind, bool diag_grp, bool diag_include,
              const double& lambda, const double& alpha, const double& gamma,
              double rho, const double tau_incr, const double tau_decr, const double nu,
              const double tol_abs, const double tol_rel, const int maxiter) {
  if (penalty == "lasso") {
    return ADMMsggl(S, group_idx, diag_ind, diag_grp, diag_include, lambda, alpha,
                    rho, tau_incr, tau_decr, nu,
                    tol_abs, tol_rel, maxiter);
  } else {
    return ADMMsggn(S, group_idx, penalty, diag_ind, diag_grp, diag_include, lambda, alpha, gamma,
                    rho, tau_incr, tau_decr, nu,
                    tol_abs, tol_rel, maxiter);
  }
}





// ADMMsggm_CV
// [[Rcpp::export]]
List ADMMsggm_CV(const arma::mat& X, const List& group_idx, std::string penalty,
                 bool diag_ind, bool diag_grp, bool diag_include,
                 const arma::colvec& lambdas, const arma::colvec& alphas, const double& gamma,
                 double rho, const double tau_incr, const double tau_decr, const double nu,
                 const double tol_abs, const double tol_rel, const int maxiter,
                 int kfold) {

  int n = X.n_rows;
  int nlambda = lambdas.n_elem;

  arma::mat CV_loss(nlambda, kfold, arma::fill::zeros);
  arma::uvec index = index_kfold(n, kfold);

  for (int k = 0; k < kfold; ++k) {

    // indices for test and training sets; test and training sets
    arma::uvec idx_test = arma::find(index == (unsigned)k);
    arma::mat X_test  = X.rows(idx_test);
    int n_test = X_test.n_rows;
    arma::uvec idx_train = arma::find(index != (unsigned)k);
    arma::mat X_train = X.rows(idx_train);

    // empirical matrix
    arma::mat S_test = arma::cov(X_test, 1);
    arma::mat S_train = arma::cov(X_train, 1);

    for (int i = 0; i < nlambda; ++i) {

      double lambda = lambdas[i];
      double alpha = alphas[i];

      List result = ADMMsggm(S_train, group_idx, penalty,
                             diag_ind, diag_grp, diag_include,
                             lambda, alpha, gamma,
                             rho, tau_incr, tau_decr, nu,
                             tol_abs, tol_rel, maxiter);

      // negative Gaussian log-likelihood
      arma::mat hatOmega = result["hatOmega"];
      double sign, logdet;
      arma::log_det(logdet, sign, hatOmega);
      if (sign <= 0 || !std::isfinite(logdet)) {
        CV_loss(i,k) = arma::datum::inf;
      } else {
        CV_loss(i,k) = (n_test / 2.0) * (arma::accu(hatOmega % S_test) - logdet);
      }

    }
  }

  // the mean of the k-fold loss for each parameter grid value
  arma::mat loss_avg = arma::mean(CV_loss, 1);
  arma::uword opt_idx = loss_avg.index_min();
  double loss_min = loss_avg(opt_idx);
  double lambda_opt = lambdas[opt_idx];
  double alpha_opt = alphas[opt_idx];

  return List::create(Named("lambda.opt") = lambda_opt,
                      Named("alpha.opt") = alpha_opt,
                      Named("CV.loss") = CV_loss,
                      Named("loss.min") = loss_min);

}


// ADMMsggm_IC
// [[Rcpp::export]]
List ADMMsggm_IC(const arma::mat& S, const List& group_idx, std::string penalty,
                 bool diag_ind, bool diag_grp, bool diag_include,
                 const arma::colvec& lambdas, const arma::colvec& alphas, const double& gamma,
                 double rho, const double tau_incr, const double tau_decr, const double nu,
                 const double tol_abs, const double tol_rel, const int maxiter,
                 std::string crit, int n, double ebic_tuning) {

  int nlambda = lambdas.n_elem;
  arma::vec IC_score(nlambda, arma::fill::zeros);
  List result_list(nlambda);

  for (int i = 0; i < nlambda; ++i) {

    double lambda = lambdas[i];
    double alpha = alphas[i];
    List result = ADMMsggm(S, group_idx, penalty,
                           diag_ind, diag_grp, diag_include,
                           lambda, alpha, gamma,
                           rho, tau_incr, tau_decr, nu,
                           tol_abs, tol_rel, maxiter);
    result_list[i] = result;

    arma::mat hatOmega = result["hatOmega"];
    IC_score[i] = criterion(hatOmega, S, n, crit, ebic_tuning);

  }

  // the mean of the score for each parameter grid value
  arma::uword opt_idx = IC_score.index_min();
  double score_min = IC_score(opt_idx);
  double lambda_opt = lambdas[opt_idx];
  double alpha_opt = alphas[opt_idx];

  return List::create(Named("lambda.opt") = lambda_opt,
                      Named("alpha.opt") = alpha_opt,
                      Named("IC.score") = IC_score,
                      Named("score.min") = score_min,
                      Named("result") = result_list[opt_idx]);

}

