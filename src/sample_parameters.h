#ifndef SAMPLE_PARAMETERS_H
#define SAMPLE_PARAMETERS_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

void sample_sig2(arma::vec& sig2_out, arma::vec& Em_str, const double a_i, const double b_i, const double T);

void res_protector(double& x);

double sample_lambda2(arma::mat& V, const double& tau, const double d_lambda, const double e_lambda, const int d,
                      const double prodlambda);

void sample_theta(arma::mat& tau2, arma::mat& coef, arma::mat& prior, const double& lambda2, const double& tau, const int r, const int c, bool Hmat);

double tau_post(double& tau, double& lambda, vec& theta, const int& d, double rat);

void sample_tau(double& tau, double& lambda, vec& theta, double& tuning, double& accept, int nburn, int irep);

double draw_bernoulli(double p);

void get_shrink(mat& V, const mat& V_prop, double& shrink, const double shrink_prop, mat& A_draw, mat& A_prior, int& accept, double& scale, int irep, int burnin);

#endif