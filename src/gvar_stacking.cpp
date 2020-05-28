// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]
#include <RcppArmadillo.h>
#include <stdlib.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "helper.h"

using namespace Rcpp;
using namespace arma;

//' @name gvar_stacking
//' @noRd
//[[Rcpp::export]]
List gvar_stacking(const SEXP xglobal_in, const SEXP plag_in, const SEXP globalpost_in, const SEXP saves_in, const SEXP thin_in,
                   const SEXP trend_in, const SEXP eigen_in) {
  //----------------------------------------------------------------------------------------------------------------------
  // GET INPUTS
  //----------------------------------------------------------------------------------------------------------------------
  NumericMatrix xglobalr(xglobal_in);
  int bigT = xglobalr.nrow(), bigK = xglobalr.ncol();
  mat xglobal(xglobalr.begin(), bigT, bigK, false);
  List globalpost(clone(globalpost_in));
  const int N = globalpost.size();
  
  const int p      = as<int>(plag_in);
  const int saves  = as<int>(saves_in);
  const int thin   = as<int>(thin_in);
  const bool trend = as<bool>(trend_in);
  const bool eigen = as<bool>(eigen_in);
  
  const int thinsaves = saves/thin;
  vec F_eigen(thinsaves, fill::zeros);
  
  int number_determinants = 1; // cons
  if(trend){number_determinants += 1;}
  
  cube A_large(thinsaves,bigK,bigK*p+number_determinants);
  cube S_large(thinsaves,bigK,bigK);
  cube Ginv_large(thinsaves,bigK,bigK);
  cube F_large(thinsaves,bigK,bigK*p);
  
  vec a1; vec b1;
  //---------------------------------------------------------------------------------------------
  vec prog_rep_points = round(linspace(0, thinsaves, 50));
  bool display_progress = true;
  Progress prog(50, display_progress);
  for(int irep = 0; irep < thinsaves; irep++){
    // patient 0
    List VAR = globalpost[0];
    mat Y = VAR["Y"];
    mat W = VAR["W"];
    unsigned int M = Y.n_cols;
    
    List store   = VAR["store"];
    List Phi     = store["Phistore"];
    List Lambda  = store["Lambdastore"];
    cube Lambda0 = store["Lambda0store"];
    cube Sigma   = store["SIGMAmed_store"];
    mat a0store  = store["a0store"];
    
    mat Lambda0irep = Lambda0.row(irep);
    if(Lambda0irep.n_cols!=M){
      Lambda0irep = Lambda0irep.t();
    }
    vec a0 = a0store.row(irep).t();
    mat A0 = join_rows(eye(M,M),-Lambda0irep.t());
    mat G  = A0*W;
    mat S  = Sigma.row(irep);
    if(trend){
      mat a1store = store["a1store"];
      a1 = a1store.row(irep).t();
    }
    List H(p);
    for(int pp=0; pp<p; pp++){
      cube Lambda_p = Lambda[pp]; mat Lambdairep = Lambda_p.row(irep);
      cube Phi_p  = Phi[pp]; mat Phiirep = Phi_p.row(irep);
      if(Lambdairep.n_cols!=M){
        Lambdairep = Lambdairep.t();
      }
      if(Phiirep.n_cols!=M){
        Phiirep = Phiirep.t();
      }
      mat H0 = join_rows(Phiirep.t(),Lambdairep.t())*W;
      H[pp]=H0;
    }
    // all others -- start at 1
    for(int cc = 1; cc < N; cc++){
      List VAR = globalpost[cc];
      mat Y = VAR["Y"];
      mat W = VAR["W"];
      unsigned int M = Y.n_cols;
      
      List store   = VAR["store"];
      List Phi     = store["Phistore"];
      List Lambda  = store["Lambdastore"];
      cube Lambda0 = store["Lambda0store"];
      cube Sigma   = store["SIGMAmed_store"];
      mat a0store  = store["a0store"];
      
      mat Lambda0irep = Lambda0.row(irep);
      if(Lambda0irep.n_cols!=M){
        Lambda0irep = Lambda0irep.t();
      }
      vec a01 = a0store.row(irep).t();
      mat A1  = join_rows(eye(M,M),-Lambda0irep.t());
      mat G1  = A1*W;
      mat S1  = Sigma.row(irep);
      // join
      G  = join_cols(G,G1);
      a0 = join_cols(a0,a01);
      if(trend){
        mat a1store = store["a1store"];
        vec a11 = a1store.row(irep).t();
        a1 = join_cols(a1,a11);
      }
      mat S0(S.n_cols,M,fill::zeros);
      S  = join_rows(S,S0);
      S1 = join_rows(S0.t(),S1);
      S  = join_cols(S,S1);
      for(int pp=0; pp<p; pp++){
        cube Lambda_p = Lambda[pp]; mat Lambdairep = Lambda_p.row(irep);
        cube Phi_p  = Phi[pp]; mat Phiirep = Phi_p.row(irep);
        if(Lambdairep.n_cols!=M){
          Lambdairep = Lambdairep.t();
        }
        if(Phiirep.n_cols!=M){
          Phiirep = Phiirep.t();
        }
        mat H1 = join_rows(Phiirep.t(),Lambdairep.t())*W;
        mat H0 = H[pp];
        mat H2 = join_cols(H0,H1);
        H[pp] = H2;
      }
    }
    mat Ginv = G.i();
    vec b0   = Ginv*a0;
    if(trend){b1 = Ginv*a1;}
    mat F; mat A;
    for(int pp=0; pp<p; pp++){
      mat temp = H[pp];
      F = join_rows(F,Ginv*temp); A = join_rows(A,Ginv*temp);
    }
    A = join_rows(A,b0);
    if(trend){A = join_rows(A,b1);}
    // save
    F_large.row(irep) = F;
    A_large.row(irep) = A;
    S_large.row(irep) = S;
    Ginv_large.row(irep) = Ginv;
    // compute eigenvalues
    if(eigen){
      mat MM(bigK*p,bigK*p, fill::zeros); MM.submat(0,0,bigK-1,bigK*p-1) = F;
      if(p>1) MM.submat(bigK,0,bigK*p-1,bigK*p-bigK-1).eye();
      cx_vec eigval; cx_mat eigvec; eig_gen(eigval, eigvec, MM);
      F_eigen(irep) = abs(real(eigval)).max();
    }
    
    // Increment progress bar
    if (any(prog_rep_points == irep)){
      prog.increment();
    }
    // check user interruption
    if(irep % 200 == 0)
      Rcpp::checkUserInterrupt();
  }
  //---------------------------------------------------------------------------------------------
  return List::create(Named("F_large")=F_large,
                      Named("S_large")=S_large,
                      Named("Ginv_large")=Ginv_large,
                      Named("A_large")=A_large,
                      Named("F_eigen")=F_eigen);
}

//' @name globalLik
//' @noRd
//[[Rcpp::export]]
List globalLik(const SEXP Y_in, const SEXP X_in, const arma::cube A_in, const arma::cube S_in, const arma::cube Ginv_in, const SEXP thinsaves_in) {
  //----------------------------------------------------------------------------------------------------------------------
  // GET INPUTS
  //----------------------------------------------------------------------------------------------------------------------
  NumericMatrix Yr(Y_in);
  NumericMatrix Xr(X_in);
  int bigT = Yr.nrow(), bigK = Yr.ncol(), bigKK = Xr.ncol();
  mat Y(Yr.begin(), bigT, bigK, false);
  mat X(Xr.begin(), bigT, bigKK, false);
  
  const int thinsaves  = as<int>(thinsaves_in);
  vec globalLik(thinsaves, fill::zeros);
  //----------------------------------------------------------------------------------------------------------------------
  // Evaluate density
  //----------------------------------------------------------------------------------------------------------------------
  for(int irep = 0; irep < thinsaves; irep++){
    mat A       = A_in.row(irep);
    mat S       = S_in.row(irep);
    mat Ginv    = Ginv_in.row(irep);
    mat Sig     = Ginv*S*Ginv.t();
    mat mean    = X*A.t();
    vec logLik  = dmvnrm_arma_fast(Y, mean, Sig, true);
    globalLik.row(irep) = sum(logLik);
  }
  return List::create(Named("globalLik")=globalLik);
}
