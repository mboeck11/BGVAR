// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include "do_rgig1.h"
using namespace Rcpp;
using namespace arma;

void sample_arcoefs(arma::mat& A_out, arma::mat& H_out, arma::mat& Em_out, arma::mat& Em_str_out,
                    arma::mat& Y, arma::mat& X, arma::mat& Sv, 
                    const arma::mat aprior, const arma::mat Vprior, const arma::mat hprior, const arma::mat Hprior) {
  // get dimensions
  int M = Y.n_cols;
  int k = X.n_cols;
  
  // Import Rs chol function
  Environment base = Environment("package:base");
  Function Rchol = base["chol"];
  
  // build Vinvprior
  cube Vinvprior(k, k, M);
  for(int mm = 0; mm < M; mm++){
    Vinvprior.slice(mm) = diagmat(1/Vprior.col(mm));
  }
  
  // estimate equation-by-equation
  for(int mm = 0; mm < M; mm++){
    if(mm == 0){
      mat S_m = exp(-0.5*Sv.col(mm));
      mat Y_m = Y.col(mm) % S_m;
      mat X_m = X % repmat(S_m,1,k);
      mat Vinv_m = Vinvprior.slice(mm);
      colvec a_m = aprior.col(mm);
      
      mat V_p = (X_m.t() * X_m + Vinv_m).i();
      mat A_p = V_p * (X_m.t() * Y_m + Vinv_m * a_m);
      
      colvec rand_normal(k);
      for(int i=0; i<k; i++){
        rand_normal(i) = R::rnorm(0,1);
      }
      mat V_p_chol_lower;
      bool chol_success = chol(V_p_chol_lower, V_p, "lower");
      // Fall back on Rs chol if armadillo fails (it suppports pivoting)
      if(chol_success == false){
        NumericMatrix tmp = Rchol(V_p, true, false, -1);
        int d = V_p.n_cols;
        mat cholV_tmp = mat(tmp.begin(), d, d, false);
        uvec piv = sort_index(as<vec>(tmp.attr("pivot")));
        V_p_chol_lower = cholV_tmp.cols(piv);
        V_p_chol_lower = V_p_chol_lower.t();
      }
      colvec A_m = A_p + V_p_chol_lower*rand_normal;
      
      A_out.col(mm) = A_m;
      Em_out.col(mm) = Y.col(mm) - X * A_m;
      Em_str_out.col(mm) = Y.col(mm) - X * A_m;
    }else{
      mat S_m = exp(-0.5*Sv.col(mm));
      mat Y_m = Y.col(mm) % S_m;
      mat X_m = join_rows(X,Em_out.cols(0,mm-1)) % repmat(S_m,1,k+mm);
      
      mat Vinv_m(k+mm, k+mm, fill::zeros);
      Vinv_m.submat(0,0,k-1,k-1) = Vinvprior.slice(mm);
      for(int i=k;i<(k+mm);i++){
        Vinv_m(i,i) = 1/Hprior(mm,i-k);
      }
      colvec a_m = join_cols(aprior.col(mm),hprior.submat(mm,0,mm,mm-1).t());
      
      mat V_p = (X_m.t() * X_m + Vinv_m).i();
      mat A_p = V_p * (X_m.t() * Y_m + Vinv_m * a_m);
      
      colvec rand_normal(k+mm);
      for(int i=0; i<(k+mm); i++){
        rand_normal(i) = R::rnorm(0,1);
      }
      mat V_p_chol_lower;
      bool chol_success = chol(V_p_chol_lower, V_p,"lower");
      // Fall back on Rs chol if armadillo fails (it suppports pivoting)
      if(chol_success == false){
        NumericMatrix tmp = Rchol(V_p, true, false, -1);
        int d = V_p.n_cols;
        mat cholV_tmp = mat(tmp.begin(), d, d, false);
        uvec piv = sort_index(as<vec>(tmp.attr("pivot")));
        V_p_chol_lower = cholV_tmp.cols(piv);
        V_p_chol_lower = V_p_chol_lower.t();
      }
      colvec A_m = A_p + V_p_chol_lower*rand_normal;
      
      A_out.col(mm) = A_m.rows(0,k-1);
      H_out.submat(mm,0,mm,mm-1) = A_m.rows(k,k+mm-1).t();
      Em_out.col(mm) = Y.col(mm) - X * A_m.rows(0,k-1);
      Em_str_out.col(mm) = Y.col(mm) - join_rows(X,Em_out.cols(0,mm-1)) * A_m;
    }
  }
}


void sample_sig2(arma::vec& sig2_out, arma::vec& Em_str, const double a_i, const double b_i, const double T){
  double a_full = a_i + T/2;
  double b_full = b_i + 0.5 * as_scalar(Em_str.t() * Em_str);
  double sig2 = 1/R::rgamma(a_full, 1/b_full);
  sig2_out.fill(sig2);
}

void res_protector(double& x){
  if (std::abs(x) < DBL_MIN * std::pow(10, 10)){
    double sign = std::copysign(1, x);
    x = DBL_MIN * std::pow(10, 10) * sign;
  }
}

double sample_lambda2(arma::mat& V, const double& tau, const double d_lambda, const double e_lambda, const int d,
                      const double prodlambda){
  double dl = d_lambda + tau*d;
  double el = e_lambda + 0.5*tau*accu(V)*prodlambda;
  double lambda2 = R::rgamma(dl, 1/el);
  res_protector(lambda2);
  return lambda2;
}

void sample_theta(arma::mat& tau2, arma::mat& coef, arma::mat& prior, const double& lambda2, const double& tau, const int r, int c, bool Hmat){
  int k;
  if(Hmat) k=1; else k=0;
  for(int i=k;i < r; i++){
    if(Hmat)  c = i;
    for(int j=0; j < c; j++){
      double lambda = tau - 0.5;
      double psi = tau * lambda2;
      double chi = std::pow(coef(i,j)-prior(i,j),2);
      
      double res = do_rgig1(lambda, chi, psi);
      if(res<1e-7) res = 1e-7;
      
      res_protector(res);
      
      tau2(i,j) = res;
    }
  }
}

double tau_post(double& tau, double& lambda, arma::vec& theta, double rat){
  double priorval = R::dexp(tau, rat, true);
  int d = theta.n_elem;
  double postval = 0;
  for(int dd=0; dd<d; dd++){
    postval = postval + R::dgamma(theta(dd), tau, 1/(tau*lambda/2), true);
  }
  double logpost = priorval + postval;
  return logpost;
}

void sample_tau(double& tau, double& lambda, arma::vec& theta, double& tuning, double& accept, int nburn, int irep){
  double old_value = tau;
  double proposal = exp(R::rnorm(0,tuning))*old_value;
  double unif = R::runif(0,1);
  
  double post_tau_prop = tau_post(proposal, lambda, theta, 1);
  double post_tau_curr = tau_post(old_value, lambda, theta, 1);
  
  double diff = post_tau_prop - post_tau_curr + std::log(proposal) - std::log(old_value);
  if(diff > log(unif)){
    tau = proposal;
    accept += 1;
  }
  if(irep < 0.5*nburn){
    if(accept/irep > 0.30){tuning = tuning*1.01;}
    if(accept/irep < 0.15){tuning = tuning*0.99;}
  }
}

double draw_bernoulli(double p){
  double unif = R::runif(0,1);
  double ret = 1;
  if(unif < p) {ret = 0;}
  return ret;
}

void get_shrink(mat& V, const mat& V_prop, double& shrink, const double shrink_prop, mat& A_draw, mat& A_prior, int& accept, double& scale, int irep, int burnin){
  int k = V.n_rows; int M = V.n_cols;
  mat postval_prop(k,M); 
  mat postval_old(k,M);
  
  // likelihood of each coefficient
  for(int i=0; i<k; i++){
    for(int j=0; j<M; j++){
      postval_prop(i,j) = log(1/sqrt(V_prop(i,j)) * exp(-0.5*pow(A_draw(i,j)-A_prior(i,j),2)/V_prop(i,j)));
      postval_old(i,j)  = log(1/sqrt(V(i,j)) * exp(-0.5*pow(A_draw(i,j)-A_prior(i,j),2)/V(i,j)));
    }
  }
  // total likelihood 
  double post_prop = accu(postval_prop) + R::dgamma(shrink_prop,0.01,1/0.01,true);  // add prior - shape scale parameterization!!!!
  double post_old = accu(postval_old) + R::dgamma(shrink,0.01,1/0.01,true); // add prior - shape scale parameterization!!!!
  if((post_prop-post_old) > log(R::runif(0,1))){
    shrink = shrink_prop;
    V = V_prop;
    accept += 1;
  }
  if((irep+1) < 0.5*burnin){
    if(accept/(irep+1) > 0.30){scale = scale*1.01;}
    if(accept/(irep+1) < 0.15){scale = scale*0.99;}
  }  
}
 

 
 
 
 
 
 
 
 
 
