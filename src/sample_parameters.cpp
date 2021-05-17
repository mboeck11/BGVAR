// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include "do_rgig1.h"
using namespace Rcpp;
using namespace arma;


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
    if(Hmat) c = i;
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
 

 
 
 
 
 
 
 
 
 
