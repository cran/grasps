#pragma once


static const double eps = 1e-10;

arma::mat soft_matrix(const arma::mat& M, const arma::mat& Thres);

arma::uvec index_kfold(const int n, const int kfold);

double criterion(const arma::mat& hatOmega, const arma::mat& S, int n,
                 const std::string& crit, double ebic_tuning);

arma::mat lambda_individual(const arma::mat& Omega, std::string penalty,
                            double lambda, double gamma);

double lambda_group(const arma::mat& Omega, std::string penalty,
                    double lambda, double gamma);

arma::mat update_Omega(const arma::mat& M, double rho);

arma::mat update_Zblock(const arma::mat& M, const std::string& penalty,
                        const Rcpp::List& group_idx,
                        const double& lambda, const double& gamma,
                        const arma::mat& initial, const double& alpha,
                        const bool diag_grp, const bool diag_include, const double rho);

