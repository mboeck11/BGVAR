// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;

arma::mat mlag(const mat& X, const int lag, const int bigT, const int bigM) {
  // get dimensions
  //int bigT = X.n_rows;
  //int bigM = X.n_cols;
  
  mat Xlag(bigT,bigM*lag,fill::zeros);
  
  for(int pp=0; pp<lag; pp++){
    int colstart = pp*bigM;
    int colend   = (pp+1)*bigM-1;
    Xlag.submat(lag,colstart,bigT-1,colend) =  X.submat(lag-pp-1,0,bigT-pp-2,bigM-1);
  }
  return Xlag;
}

arma::vec seq(const int first, const int last){
  NumericVector y(abs(last-first)+1);
  if(first<last){
    std::iota(y.begin(), y.end(), first);
  }else{
    std::iota(y.begin(), y.end(), last);
    std::reverse(y.begin(), y.end());
  }
  vec ret(y.begin(), abs(last-first)+1, false);
  return ret;
}

void gen_compMat(mat& Cm, const mat& A, const int M, const int p) {
  int Mp = M*p;
  Cm.submat(0,0,M-1,Mp-1) = A.t();
  if(p>1) Cm.submat(M,0,Mp-1,Mp-M-1).eye();
}

double get_ar(mat& Yraw, int p){
  int Traw = Yraw.n_rows;
  mat Xraw = mlag(Yraw,p,Traw,1);
  mat X = Xraw.submat(p,0,Traw-1,0);
  mat Y = Yraw.submat(p,0,Traw-1,0);
  
  int T = Y.n_rows;
  
  vec trendvec(T); for(int tt=0; tt<T; tt++) {trendvec(tt)=tt+1;} 
  X = join_rows(X,trendvec);
  
  vec b_OLS = inv(X.t() * X) * X.t() * Y;
  vec e_OLS = Y - X * b_OLS;
  double s_OLS = dot(e_OLS.t(), e_OLS) / (T-p-1);
  
  return s_OLS;
}

void get_Vminnesota(mat& V, vec& sigmas, double shrink1, double shrink2, double shrink3, double shrink4, bool cons, int Mstar, 
                    int plag, int plagstar, bool trend){
  int k = V.n_rows; int M = V.n_cols;
  
  // endogenous part
  for(int i=0; i < M; i++){ // equation - column
    for(int pp=1; pp <= plag; pp++){
      for(int j=0; j < M; j++){ // variable - row
        if(i==j){ // own lag
          V(j+M*(pp-1),i) = pow(shrink1/static_cast<double>(pp),2.0);
        }else{
          V(j+M*(pp-1),i) = pow((shrink1 * shrink2)/(static_cast<double>(pp)),2.0) * (sigmas(i)/sigmas(j));
        }
      }
    }
  }
  // exogenous part
  if(Mstar > 0){
    for(int i=0; i < M; i++){ // equation - column
      for(int pp=0; pp <= plagstar; pp++){
        for(int j=0; j < Mstar; j++){
          V(M*plag+pp*Mstar+j,i) = pow((shrink1 * shrink3)/static_cast<double>(pp+1),2.0) * sigmas(i)/(sigmas(M+j));
        }
      }
    }
  }
  // deterministics
  if(cons){
    for(int i=0; i < M; i++){
      V(k-1,i) = shrink4*sigmas(i);
      if(trend) V(k-2,i) = shrink4*sigmas(i);
    }
  }
}

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

//' @name dmvnrm_arma_fast
//' @noRd
//[[Rcpp::export]]
arma::vec dmvnrm_arma_fast(const arma::mat& x,  
                           const arma::mat& mean,  
                           const arma::mat& sigma, 
                           bool const logd = false) { 
  // Import Rs chol function
  Environment base = Environment("package:base");
  Function Rchol = base["chol"];
  
  double const log2pi = std::log(2.0 * M_PI);
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  
  // fall back on pivoting if Cholesky do not work
  mat Sig_chol;
  bool chol_success = chol(Sig_chol, sigma);
  // Fall back on Rs chol if armadillo fails (it suppports pivoting)
  if(chol_success == false){
    NumericMatrix tmp = Rchol(sigma, true);
    int d = sigma.n_cols;
    mat cholV_tmp = mat(tmp.begin(), d, d, false);
    uvec piv = sort_index(as<vec>(tmp.attr("pivot")));
    Sig_chol = cholV_tmp.cols(piv);
  }
  
  arma::mat const rooti = arma::inv(trimatu(Sig_chol)); // arma::chol(sigma)
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean.row(i));
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd){
    return out;
  }else{
    return exp(out);
  }
}

arma::vec dmvnrm_arma_old(mat& x,  
                          mat& mean,  
                          mat& sigma, 
                          bool logd = false) { 
  double const log2pi = std::log(2.0 * M_PI);
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(double)xdim/2.0 * log2pi;
  
  for (uword i = 0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean.row(i) ) ;    
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
  }  
  
  if (logd)
    return out;
  return exp(out);
}

arma::mat robust_chol(const arma::mat& V){
  // Note that this function returns the lower triangular form,
  // i.e. cholV * cholV.t() = V;
  
  arma::mat cholV;
  bool chol_success = arma::chol(cholV, V);
  
  if (chol_success == false) {
    int max_tries = 1000;
    int num_tries = 1;
    double jitter = 1e-12 * arma::mean(V.diag());
    // Jitter main diagonal to (potentially) solve numerical issues
    while ((chol_success == false) & (num_tries <= max_tries) & !std::isinf(jitter)){
      chol_success = arma::chol(cholV, V + jitter * arma::eye(arma::size(V)));
      
      jitter *= 1.1;
      num_tries += 1;
    }
  }
  
  if (chol_success == true) {
    return cholV.t();
  } else {
    return arma::mat(arma::size(V), arma::fill::none);
  }
}
