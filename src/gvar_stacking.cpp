// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdlib.h>
#include "helper.h"

using namespace Rcpp;
using namespace arma;

//' @name gvar_stacking
//' @noRd
//[[Rcpp::export]]
List gvar_stacking(const arma::mat xglobal, const int plag, const Rcpp::List globalpost, const int draws, const int thin,
                   const bool trend, const bool eigen, const bool verbose) {
  //----------------------------------------------------------------------------------------------------------------------
  // GET INPUTS
  //----------------------------------------------------------------------------------------------------------------------
  const int bigK = xglobal.n_cols;
  const int N = globalpost.size();
  
  const int thindraws = draws/thin;
  vec F_eigen(thindraws, fill::zeros);
  
  int number_determinants = 1; // cons
  if(trend){number_determinants += 1;}
  
  arma::cube A_large(bigK,bigK*plag+number_determinants,thindraws);
  arma::cube S_large(bigK,bigK,thindraws);
  arma::cube Ginv_large(bigK,bigK,thindraws);
  arma::cube F_large(bigK,bigK*plag,thindraws);
  
  arma::vec a1, b1;
  //---------------------------------------------------------------------------------------------
  //vec prog_rep_points = round(linspace(0, thindraws, 50));
  //bool display_progress = true;
  //Progress prog(50, verbose);
  for(int irep = 0; irep < thindraws; irep++){
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
    
    mat Lambda0irep = Lambda0.slice(irep);
    arma::colvec a0 = a0store.col(irep);
    arma::mat A0 = join_rows(eye(M,M),-Lambda0irep.t());
    arma::mat G  = A0*W;
    arma::mat S  = Sigma.slice(irep);
    if(trend){
      mat a1store = store["a1store"];
      a1 = a1store.col(irep);
    }
    List H(plag);
    for(int pp=0; pp < plag; pp++){
      arma::cube Lambda_p = Lambda[pp]; 
      arma::cube Phi_p = Phi[pp]; 
      arma::mat Lambdairep = Lambda_p.slice(irep);
      arma::mat Phiirep = Phi_p.slice(irep);
      
      mat H0 = join_rows(Phiirep.t(),Lambdairep.t())*W;
      H[pp]=H0;
    }
    // all others -- start at 1
    for(int cc = 1; cc < N; cc++){
      List VAR = globalpost[cc];
      arma::mat Y = VAR["Y"];
      arma::mat W = VAR["W"];
      unsigned int M = Y.n_cols;
      
      Rcpp::List store   = VAR["store"];
      Rcpp::List Phi     = store["Phistore"];
      Rcpp::List Lambda  = store["Lambdastore"];
      arma::cube Lambda0 = store["Lambda0store"];
      arma::cube Sigma   = store["SIGMAmed_store"];
      arma::mat a0store  = store["a0store"];
      
      arma::mat Lambda0irep = Lambda0.slice(irep);
      arma::colvec a01 = a0store.col(irep);
      arma::mat A1  = join_rows(eye(M,M),-Lambda0irep.t());
      arma::mat G1  = A1*W;
      arma::mat S1  = Sigma.slice(irep);
      // join
      G  = join_cols(G,G1);
      a0 = join_cols(a0,a01);
      if(trend){
        mat a1store = store["a1store"];
        vec a11 = a1store.col(irep);
        a1 = join_cols(a1,a11);
      }
      mat S0(S.n_cols,M,fill::zeros);
      S  = join_rows(S,S0);
      S1 = join_rows(S0.t(),S1);
      S  = join_cols(S,S1);
      for(int pp=0; pp < plag; pp++){
        arma::cube Lambda_p = Lambda[pp]; 
        arma::cube Phi_p = Phi[pp]; 
        arma::mat Lambdairep = Lambda_p.slice(irep);
        arma::mat Phiirep = Phi_p.slice(irep);
       
        mat H1 = join_rows(Phiirep.t(),Lambdairep.t())*W;
        mat H0 = H[pp];
        mat H2 = join_cols(H0,H1);
        H[pp] = H2;
      }
    }
    
    mat Ginv = G.i();
    vec b0   = Ginv*a0;
    if(trend){b1 = Ginv*a1;}
    arma::mat F, A;
    for(int pp=0; pp < plag; pp++){
      arma::mat temp = H[pp];
      F = join_rows(F,Ginv*temp); A = join_rows(A,Ginv*temp);
    }
    A = join_rows(A,b0);
    if(trend){A = join_rows(A,b1);}
    
    
    // save
    F_large.slice(irep) = F;
    A_large.slice(irep) = A;
    S_large.slice(irep) = S;
    Ginv_large.slice(irep) = Ginv;
  
    // compute eigenvalues
    if(eigen){
      arma::mat MM(bigK*plag, bigK*plag, fill::zeros); MM.submat(0,0,bigK-1,bigK*plag-1) = F;
      if(plag>1) MM.submat(bigK,0,bigK*plag-1,bigK*plag-bigK-1).eye();
      cx_vec eigval; cx_mat eigvec; eig_gen(eigval, eigvec, MM);
      F_eigen(irep) = abs(real(eigval)).max();
    }
    
    //if(verbose){
      // Increment progress bar
      //if (any(prog_rep_points == irep)){
      //  prog.increment();
      //}
    //}
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
List globalLik(const SEXP Y_in, const SEXP X_in, const arma::cube A_in, const arma::cube S_in, const arma::cube Ginv_in, const SEXP thindraws_in){
  //----------------------------------------------------------------------------------------------------------------------
  // GET INPUTS
  //----------------------------------------------------------------------------------------------------------------------
  NumericMatrix Yr(Y_in);
  NumericMatrix Xr(X_in);
  int bigT = Yr.nrow(), bigK = Yr.ncol(), bigKK = Xr.ncol();
  mat Y(Yr.begin(), bigT, bigK, false);
  mat X(Xr.begin(), bigT, bigKK, false);
  
  const int thindraws  = as<int>(thindraws_in);
  vec globalLik(thindraws, fill::zeros);
  //----------------------------------------------------------------------------------------------------------------------
  // Evaluate density
  //----------------------------------------------------------------------------------------------------------------------
  for(int irep = 0; irep < thindraws; irep++){
    mat A       = A_in.slice(irep);
    mat S       = S_in.slice(irep);
    mat Ginv    = Ginv_in.slice(irep);
    mat Sig     = Ginv*S*Ginv.t();
    mat mean    = X*A.t();
    vec logLik  = dmvnrm_arma_fast(Y, mean, Sig, true);
    globalLik.row(irep) = sum(logLik);
  }
  return List::create(Named("globalLik")=globalLik);
}
