#ifndef HELPER_H
#define HELPER_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

arma::mat mlag(const mat& X, const int lag, const int bigT, const int bigM);

arma::vec seq(const int first, const int last);

void gen_compMat(mat& Cm, const mat& A, const int M, const int p);

double get_ar(mat& Y, int p);

void get_Vminnesota(mat& V, vec& sigmas, double shrink1, double shrink2, double shrink3, double shrink4, bool cons, int Mstar, int p, bool trend);

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat);

arma::vec dmvnrm_arma_fast(const arma::mat& x, const arma::mat& mean, const arma::mat& sigma, bool const logd = false);

arma::vec dmvnrm_arma_old(arma::mat& x, arma::mat& mean, arma::mat& sigma, bool logd = false);

#endif