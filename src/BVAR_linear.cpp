#include <RcppArmadillo.h>
#include <stochvol.h>
#include <stdlib.h>
#include "helper.h"
#include "do_rgig1.h"

using namespace Rcpp;
using namespace arma;

double draw_bernoulli(double p){
  double unif = R::runif(0,1);
  double ret = 1;
  if(unif < p) {ret = 0;}
  return ret;
}

double tau_post(double tau, double lambda, arma::vec theta, double rat){
  double priorval = R::dexp(tau, rat, true);
  int d = theta.n_elem;
  double postval = 0;
  for(int dd=0; dd<d; dd++){
    postval = postval + R::dgamma(theta(dd), tau, 1/(tau*lambda/2), true);
  }
  double logpost = priorval + postval;
  return logpost;
}

//' @name BVAR_linear
//' @noRd
//[[Rcpp::interfaces(r, cpp)]]
//[[Rcpp::export]]
List BVAR_linear(arma::mat Yraw, 
                 arma::mat Wraw, 
                 arma::mat Exraw,
                 int plag, 
                 int draws, 
                 int burnin,
                 int thin,
                 bool cons, 
                 bool trend,
                 bool sv,
                 int prior,
                 Rcpp::List hyperparam) {
  //----------------------------------------------------------------------------------------------------------------------
  // CONSTRUCT DATA
  //----------------------------------------------------------------------------------------------------------------------
  int Traw = Yraw.n_rows;
  int M = Yraw.n_cols;
  int K = M*plag;
  int Mstar = Wraw.n_cols;
  int Kstar = Mstar*(plag+1);
  
  bool texo = false; int Mex=0;
  if(Exraw.n_elem != 1){
    texo = true; Mex = Exraw.n_cols; 
  }
  
  mat Xraw = mlag(Yraw,plag,Traw,M);
  mat X0 = Xraw.submat(plag,0,Traw-1,K-1);
  mat Y = Yraw.submat(plag,0,Traw-1,M-1);
  double T = X0.n_rows;
  mat Wall = join_rows(Wraw,mlag(Wraw,plag,Traw,Mstar));
  mat W0 = Wall.submat(plag,0,Traw-1,Kstar-1);
  mat X = join_rows(X0,W0);
  if(texo){
    mat E0 = Exraw.submat(plag,0,Traw-1,Mex-1);
    X = join_rows(X,E0);
  }
  if(cons){
    X = join_rows(X,colvec(T,fill::ones));
  }
  if(trend){
    vec trendvec(T); for(int tt=0; tt<T; tt++) {trendvec(tt)=tt+1;}
    X = join_rows(X,trendvec);
  }

  int k = X.n_cols;
  int v = (M*(M-1))/2;
  //----------------------------------------------------------------------------------------------------------------------
  // HYPERPARAMETERS
  //----------------------------------------------------------------------------------------------------------------------
  // prior == 1: SIMS
  double shrink1 = hyperparam["shrink1"];
  double shrink2 = hyperparam["shrink2"];
  double shrink3 = hyperparam["shrink3"];
  double shrink4= hyperparam["shrink4"];
  // prior == 2: SSVS
  const double tau00 = hyperparam["tau0"];
  const double tau11 = hyperparam["tau1"];
  double p_i = hyperparam["p_i"];
  const double kappa00 = hyperparam["kappa0"];
  const double kappa11 = hyperparam["kappa1"];
  double q_ij = hyperparam["q_ij"];
  // prior == 3: NG
  const double d_lambda = hyperparam["d_lambda"];
  const double e_lambda = hyperparam["e_lambda"];
  const double tau_theta = hyperparam["tau_theta"];
  const bool sample_tau_bool = hyperparam["sample_tau"];
  // misc
  const double a_1 = hyperparam["a_1"];
  const double b_1 = hyperparam["b_1"];
  const double prmean = hyperparam["prmean"];
  const double Bsigma = hyperparam["Bsigma"];
  const double a0 = hyperparam["a0"];
  const double b0 = hyperparam["b0"];
  const double bmu = hyperparam["bmu"];
  const double Bmu = hyperparam["Bmu"];
  //----------------------------------------------------------------------------------------------------------------------
  // OLS QUANTITIES
  //----------------------------------------------------------------------------------------------------------------------
  mat XtXinv = (X.t() * X).i();
  mat A_OLS = XtXinv * (X.t()) * Y;
  mat E_OLS = Y - X * A_OLS;
  mat SIGMA_OLS = ((E_OLS.t()) * E_OLS) / (T-k);
  //----------------------------------------------------------------------------------------------------------------------
  // INITIAL VALUES
  //----------------------------------------------------------------------------------------------------------------------
  // initial values and miscellaneous stuff
  mat A_draw = A_OLS;
  cube SIGMA(M,M,T);
  for(int tt=0; tt < T; tt++){
    SIGMA.slice(tt) = SIGMA_OLS;
  }
  mat Em_draw = E_OLS; mat Em_str_draw = E_OLS;
  mat L_draw(M,M); L_draw.eye();
  mat Cm(K,K, fill::zeros); gen_compMat(Cm, A_OLS.rows(0,K-1), M, plag);
  //----------------------------------------------------------------------------------------------------------------------
  // PRIORS
  //----------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------------------------
  // prior on VAR coefficients
  //---------------------------------------------------------------
  // prior mean
  mat A_prior(k,M); A_prior.fill(0); A_prior.diag() += prmean;
  
  // prior variance
  mat V_prior(k,M); V_prior.fill(10);
  
  // SIMS stuff
  double accept1 = 0, accept2 = 0, accept4 = 0;
  double scale1 = 0.43, scale2 = 0.43, scale4 = 0.43;
  vec sigmas(M+Mstar, fill::zeros);
  for(int i=0; i < M; i++){
    mat Y_ = Yraw.col(i);
    sigmas(i) = get_ar(Y_,plag);
  }
  for(int i=0; i < Mstar; i++){
    mat W_ = Wraw.col(i);
    sigmas(M+i) = get_ar(W_,plag);
  }
  if(prior==1) get_Vminnesota(V_prior, sigmas, shrink1, shrink2, shrink3, shrink4, cons, Mstar, plag, trend);
  
  // initialize stuff for MN prior
  mat V_prop1(k,M), V_prop2(k,M), V_prop4(k,M); 
  double shrink_prop1, shrink_prop2, shrink_prop4;
  double post_old1=0.0, post_prop1=0.0, post_old2=0.0, post_prop2=0.0, post_old4=0.0, post_prop4=0.0;
  for(int i=0; i<k; i++){
    for(int j=0; j<M; j++){
      post_old1 = post_old1 + R::dnorm(A_draw(i,j),A_prior(i,j),std::sqrt(V_prior(i,j)),true);
      post_old2 = post_old2 + R::dnorm(A_draw(i,j),A_prior(i,j),std::sqrt(V_prior(i,j)),true);
      post_old4 = post_old4 + R::dnorm(A_draw(i,j), A_prior(i,j), std::sqrt(V_prior(i,j)), true);
    }
  }
  post_old1 = post_old1 + R::dgamma(shrink1,0.01,1/0.01,true) + log(shrink1); // add prior - shape scale parameterization!!!! + correction term
  post_old2 = post_old2 + R::dgamma(shrink2,0.01,1/0.01,true) + log(shrink2); // add prior - shape scale parameterization!!!! + correction term
  post_old4 = post_old4 + R::dgamma(shrink4,0.01,1/0.01,true) + log(shrink4); // add prior - shape scale parameterization!!!! + correction term
  
  // SSVS stuff
  mat gamma(k,M, fill::ones);
  mat temp = kron(SIGMA_OLS,XtXinv);
  vec sigma_alpha = arma::sqrt(diagvec(kron(SIGMA_OLS,XtXinv))); // scale with OLS standard deviation
  mat tau0(k, M, fill::zeros); mat tau1(k, M, fill::zeros);
  int ii=0;
  for(int j=0; j < M; j++){
    for(int i=0; i < k; i++){
      tau0(i,j) = tau00 * sigma_alpha(ii); // tau00: small sd
      tau1(i,j) = tau11 * sigma_alpha(ii); // tau11: big sd
      ii += 1;
    }
  }
  // NG stuff
  mat lambda2_A(plag+1,2,fill::zeros);
  mat A_tau(plag+1,2); A_tau.fill(tau_theta); A_tau(0,0)=0;
  mat A_tuning(plag+1,2); A_tuning.fill(0.43);
  mat A_accept(plag+1,2, fill::zeros);
  // initialize stuff for NG prior
  mat A_con, V_con, P_con, A_end, V_end, P_end, A_exo, V_exo, P_exo;
  int r_con, c_con, d_con, r_end, c_end, d_end, r_exo, c_exo, d_exo;
  double prodlambda, dl, el, lambda, chi, psi, res;
  double unif, proposal, post_tau_prop, post_tau_curr, diff;
  //---------------------------------------------------------------
  // prior on coefficients in H matrix of VCV
  //---------------------------------------------------------------
  // prior mean
  mat l_prior(M,M,fill::zeros);

  // prior variance
  mat L_prior(M,M); L_prior.fill(kappa11); L_prior = trimatl(L_prior); L_prior.diag().zeros();
  
  // SSVS
  mat omega(M,M, fill::ones); omega = trimatl(omega); omega.diag().zeros();
  
  // NG
  double lambda2_L = 0.01;
  double L_tau = tau_theta;
  double L_tuning = 0.43;
  double L_accept = 0.0;
  mat lambda2_Lmat(plag+1,1, fill::zeros);
  mat L_taumat(plag+1,1, fill::zeros); 
  mat L_accmat(plag+1,1, fill::zeros);
  mat L_tunmat(plag+1,1, fill::zeros);
  //---------------------------------------------------------------
  // SV quantitites
  //---------------------------------------------------------------
  mat Sv_draw(T,M); Sv_draw.fill(-3);
  mat Sv_para(4,M);
  for(int mm=0; mm < M; mm++){
    Sv_para(0,mm) = -10;
    Sv_para(1,mm) = .9;
    Sv_para(2,mm) = .2;
    Sv_para(3,mm) = -10;
  }
  uvec rec(T); rec.fill(5);
  const double offset = 0;  // maybe want to change to 1e-40 or so to be on the safe side? I have got random NAs because of log(0) for real data inputs
  using stochvol::PriorSpec;
  const PriorSpec prior_spec = {  // prior specification object for the update_*_sv functions
    PriorSpec::Latent0(),  // stationary prior distribution on priorlatent0
    PriorSpec::Mu(PriorSpec::Normal(bmu, std::sqrt(Bmu))),  // normal prior on mu
    PriorSpec::Phi(PriorSpec::Beta(a0, b0)),  // stretched beta prior on phi
    PriorSpec::Sigma2(PriorSpec::Gamma(0.5, 0.5 / Bsigma))  // normal(0, Bsigma) prior on sigma
  };  // heavy-tailed, leverage, regression turned off
  using stochvol::ExpertSpec_FastSV;
  const ExpertSpec_FastSV expert {  // very expert settings for the Kastner, Fruehwirth-Schnatter (2014) sampler
    true,  // interweave
    stochvol::Parameterization::CENTERED,  // centered baseline always
    1e-8,  // B011inv,
    1e-12,  //B022inv,
    2,  // MHsteps,
    ExpertSpec_FastSV::ProposalSigma2::INDEPENDENCE,  // independece proposal for sigma
    -1,  // unused for independence prior for sigma
    ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL  // immediately reject (mu,phi,sigma) if proposed phi is outside (-1, 1)
  };
  // initialize stuff
  double a_full, b_full, sig2, mu, phi, sigma, h0;
  //---------------------------------------------------------------------------------------------------------------
  // SAMPLER MISCELLANEOUS
  //---------------------------------------------------------------------------------------------------------------
  int ntot = burnin + draws;

  // thinning parameters
  int thindraws = draws/thin;
  
  // import R's chol function
  Environment base = Environment("package:base");
  Function Rchol = base["chol"];
  //---------------------------------------------------------------------------------------------
  // STORAGES
  //---------------------------------------------------------------------------------------------
  cube A_store(thindraws,k,M);
  cube L_store(thindraws,M,M);
  cube res_store(thindraws,T,M);
  // SV
  cube Sv_store(thindraws,T,M);
  cube pars_store(thindraws,4,M);
  // SIMS
  mat shrink_store(thindraws,3, fill::zeros);
  // SSVS
  cube gamma_store(thindraws,k,M, fill::zeros);
  cube omega_store(thindraws,M,M, fill::zeros);
  // NG
  cube theta_store(thindraws,k,M, fill::zeros);
  cube lambda2_store(thindraws,plag+1,3, fill::zeros);
  cube tau_store(thindraws,plag+1,3, fill::zeros);
  //---------------------------------------------------------------------------------------------
  // MCMC LOOP
  //---------------------------------------------------------------------------------------------
  for(int irep = 0; irep < ntot; irep++){
    // Step 1: Sample coefficients
    for(int mm = 0; mm < M; mm++){ // estimate equation-by-equation
      if(mm == 0){
        mat S_m = exp(-0.5*Sv_draw.col(mm));
        mat Y_m = Y.col(mm) % S_m;
        mat X_m = X % repmat(S_m,1,k);
        mat Vinv_m = diagmat(1/V_prior.col(mm));
        colvec a_m = A_prior.col(mm);
        
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
        
        A_draw.col(mm) = A_m;
        Em_draw.col(mm) = Y.col(mm) - X * A_m;
        Em_str_draw.col(mm) = Y.col(mm) - X * A_m;
      }else{
        mat S_m = exp(-0.5*Sv_draw.col(mm));
        mat Y_m = Y.col(mm) % S_m;
        mat X_m = join_rows(X,Em_draw.cols(0,mm-1)) % repmat(S_m,1,k+mm);
        
        mat Vinv_m(k+mm, k+mm, fill::zeros);
        Vinv_m.submat(0,0,k-1,k-1) = diagmat(1/V_prior.col(mm));
        for(int i=k;i<(k+mm);i++){
          Vinv_m(i,i) = 1/L_prior(mm,i-k);
        }
        colvec a_m = join_cols(A_prior.col(mm),l_prior.submat(mm,0,mm,mm-1).t());
        
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
        
        A_draw.col(mm) = A_m.rows(0,k-1);
        L_draw.submat(mm,0,mm,mm-1) = A_m.rows(k,k+mm-1).t();
        Em_draw.col(mm) = Y.col(mm) - X * A_m.rows(0,k-1);
        Em_str_draw.col(mm) = Y.col(mm) - join_rows(X,Em_draw.cols(0,mm-1)) * A_m;
      }
    }
    //-----------------------------------------------
    // Step 2: different prior setups
    // SIMS
    if(prior == 1){
      // first shrinkage parameter (own lags)
      shrink_prop1 = exp(R::rnorm(0,scale1))*shrink1;
      get_Vminnesota(V_prop1, sigmas, shrink_prop1, shrink2, shrink3, shrink4, cons, Mstar, plag, trend);
      // likelihood of each coefficient
      for(int i=0; i<k; i++){
        for(int j=0; j<M; j++){
          post_prop1 = post_prop1 + R::dnorm(A_draw(i,j), A_prior(i,j), std::sqrt(V_prop1(i,j)), true);
        }
      }
      // total likelihood 
      post_prop1 = post_prop1 + R::dgamma(shrink_prop1,0.01,1/0.01,true) + log(shrink_prop1);  // add prior - shape scale parameterization!!!! + correction term
      if((post_prop1-post_old1) > log(R::runif(0,1))){
        shrink1 = shrink_prop1;
        V_prior = V_prop1;
        post_old1 = post_prop1;
        accept1 += 1;
      }
      
      // second shrinkage parameter (cross equations)
      shrink_prop2 = exp(R::rnorm(0,scale2))*shrink2;
      get_Vminnesota(V_prop2, sigmas, shrink1, shrink_prop2, shrink3, shrink4, cons, Mstar, plag, trend);
      // likelihood of each coefficient
      for(int i=0; i<k; i++){
        for(int j=0; j<M; j++){
          post_prop2 = post_prop2 + R::dnorm(A_draw(i,j), A_prior(i,j), std::sqrt(V_prop2(i,j)), true);
        }
      }
      // total likelihood 
      post_prop2 = post_prop2 + R::dgamma(shrink_prop2,0.01,1/0.01,true) + log(shrink_prop2);  // add prior - shape scale parameterization!!!! + correction term
      if((post_prop2-post_old2) > log(R::runif(0,1))){
        shrink2 = shrink_prop2;
        V_prior = V_prop2;
        post_old2 = post_prop2; 
        accept2 += 1;
      }
      
      // fourth shrinkage parameter (weakly exogenous)
      shrink_prop4 = exp(R::rnorm(0,scale4))*shrink4;
      get_Vminnesota(V_prop4, sigmas, shrink1, shrink2, shrink3, shrink_prop4, cons, Mstar, plag, trend);
      // likelihood of each coefficient
      for(int i=0; i<k; i++){
        for(int j=0; j<M; j++){
          post_prop4 = post_prop4 + R::dnorm(A_draw(i,j), A_prior(i,j), std::sqrt(V_prop4(i,j)), true);
        }
      }
      // total likelihood 
      post_prop4 = post_prop4 + R::dgamma(shrink_prop4,0.01,1/0.01,true) + log(shrink_prop4);  // add prior - shape scale parameterization!!!! + correction term
      if((post_prop4-post_old4) > log(R::runif(0,1))){
        shrink4 = shrink_prop4;
        V_prior = V_prop4;
        post_old4 = post_prop4;
        accept4 += 1;
      }
      
      if((irep+1) < 0.5*burnin){
        if(accept1/(irep+1) > 0.30){scale1 *= 1.01;}
        if(accept1/(irep+1) < 0.15){scale1 *= 0.99;}
        if(accept2/(irep+1) > 0.30){scale2 *= 1.01;}
        if(accept2/(irep+1) < 0.15){scale2 *= 0.99;}
        if(accept4/(irep+1) > 0.30){scale4 *= 1.01;}
        if(accept4/(irep+1) < 0.15){scale4 *= 0.99;}
      }  
    }
    // SSVS
    if(prior == 2){
      // coefficients A matrix
      for(int j=0; j < M; j++){
        for(int i=0; i < k; i++){
          double u_i1 = R::dnorm(A_draw(i,j), A_prior(i,j), tau0(i,j),false) * p_i;
          double u_i2 = R::dnorm(A_draw(i,j), A_prior(i,j), tau1(i,j),false) * (1-p_i);
          double ast  = u_i1/(u_i1+u_i2);
          if(NumericVector::is_na(ast)) ast = 0;
          gamma(i,j) = draw_bernoulli(ast);
          if(gamma(i,j)==0) {V_prior(i,j) = tau0(i,j)*tau0(i,j);}
          if(gamma(i,j)==1) {V_prior(i,j) = tau1(i,j)*tau1(i,j);}
        }
      }
      // coefficients H matrix
      for(int i=1; i < M; i++){
        for(int j=0; j < i; j++){
          double u_ij1 = R::dnorm(L_draw(i,j),l_prior(i,j),kappa00,false) * q_ij;
          double u_ij2 = R::dnorm(L_draw(i,j),l_prior(i,j),kappa11,false) * (1-q_ij);
          double hst = u_ij1/(u_ij1+u_ij2);
          if(NumericVector::is_na(hst)) hst = 1;
          omega(i,j) = draw_bernoulli(hst);
          if(omega(i,j)==0) {L_prior(i,j) = kappa00*kappa00;}
          if(omega(i,j)==1) {L_prior(i,j) = kappa11*kappa11;}
        }
      }
    }
    // NG
    if(prior == 3){
      // coefficients A matrix
      for(int pp=0; pp < (plag+1); pp++){
        if(pp==0){
          A_con = A_draw.rows(plag*M, plag*M+Mstar-1); 
          V_con = V_prior.rows(plag*M, plag*M+Mstar-1); 
          P_con = A_prior.rows(plag*M, plag*M+Mstar-1);
          r_con = A_con.n_rows; c_con = A_con.n_cols; d_con = A_con.n_elem;
          
          // sample lambda
          dl = d_lambda + A_tau(0,1)*d_con;
          el = e_lambda + 0.5*A_tau(0,1)*accu(V_con);
          lambda2_A(0,1) = R::rgamma(dl, 1/el);
          
          // sample theta
          for(int ii=0; ii < r_con; ii++){
            for(int jj=0; jj < c_con; jj++){
              lambda = A_tau(0,1) - 0.5;
              psi = A_tau(0,1) * lambda2_A(0,1);
              chi = std::pow(A_con(ii,jj)-P_con(ii,jj),2);
              
              res = do_rgig1(lambda, chi, psi);
              if(res<1e-7) res = 1e-7;
              if(res>1e+7) res = 1e+7;
              
              V_con(ii,jj) = res;
            }
          }
          V_prior.rows(plag*M, plag*M+Mstar-1) = V_con;
          
          // sample tau
          if(sample_tau_bool){
            vec theta_vec = V_con.as_col();
            prodlambda = as_scalar(lambda2_A.submat(0,1,pp,1));
            
            proposal = exp(R::rnorm(0,A_tuning(0,1)))*A_tau(0,1);
            unif = R::runif(0,1);
            
            post_tau_prop = tau_post(proposal, prodlambda, theta_vec, 1);
            post_tau_curr = tau_post(A_tau(0,1), prodlambda, theta_vec, 1);
            diff = post_tau_prop - post_tau_curr + std::log(proposal) - std::log(A_tau(0,1));
            if(diff > log(unif)){
              A_tau(0,1) = proposal;
              A_accept(0,1) += 1;
            }
            if(irep < 0.5*burnin){
              if(A_accept(0,1)/irep > 0.30){A_tuning(0,1) = A_tuning(0,1)*1.01;}
              if(A_accept(0,1)/irep < 0.15){A_tuning(0,1) = A_tuning(0,1)*0.99;}
            }
          }
        }else{
          A_end = A_draw.rows((pp-1)*M, pp*M-1); 
          V_end = V_prior.rows((pp-1)*M, pp*M-1); 
          P_end = A_prior.rows((pp-1)*M, pp*M-1);
          r_end = A_end.n_rows; c_end = A_end.n_cols; d_end = A_end.n_elem;
          
          // sample lambda
          if(pp == 1){
            prodlambda = 1.0;
          }else{
            prodlambda = as_scalar(prod(lambda2_A.submat(1,0,pp-1,0)));
          }
          dl = d_lambda + A_tau(pp,0)*d_end;
          el = e_lambda + 0.5*A_tau(pp,0)*accu(V_end)*prodlambda;
          lambda2_A(pp,0) = R::rgamma(dl, 1/el);
          
          // sample theta
          prodlambda = as_scalar(prod(lambda2_A.submat(1,0,pp,0)));
          for(int ii=0; ii < r_end; ii++){
            for(int jj=0; jj < c_end; jj++){
              lambda = A_tau(pp,0) - 0.5;
              psi = A_tau(pp,0) * prodlambda;
              chi = std::pow(A_end(ii,jj)-P_end(ii,jj),2);
              
              res = do_rgig1(lambda, chi, psi);
              if(res<1e-7) res = 1e-7;
              if(res>1e+7) res = 1e+7;
              
              V_end(ii,jj) = res;
            }
          }
          V_prior.rows((pp-1)*M, pp*M-1) = V_end;
          
          // sample tau
          if(sample_tau_bool){
            vec theta_vec = V_end.as_col();
            prodlambda = as_scalar(lambda2_A.submat(1,0,pp,0));
            
            proposal = exp(R::rnorm(0,A_tuning(pp,0)))*A_tau(pp,0);
            unif = R::runif(0,1);
            
            post_tau_prop = tau_post(proposal, prodlambda, theta_vec, 1);
            post_tau_curr = tau_post(A_tau(pp,0), prodlambda, theta_vec, 1);
            diff = post_tau_prop - post_tau_curr + std::log(proposal) - std::log(A_tau(pp,0));
            if(diff > log(unif)){
              A_tau(pp,0) = proposal;
              A_accept(pp,0) += 1;
            }
            if(irep < 0.5*burnin){
              if(A_accept(pp,0)/irep > 0.30){A_tuning(pp,0) = A_tuning(pp,0)*1.01;}
              if(A_accept(pp,0)/irep < 0.15){A_tuning(pp,0) = A_tuning(pp,0)*0.99;}
            }
          }
          //-------------------------------------------------------------------
          // weakly exogenous
          A_exo = A_draw.rows(plag*M+pp*Mstar, plag*M+(pp+1)*Mstar-1); 
          V_exo = V_prior.rows(plag*M+pp*Mstar, plag*M+(pp+1)*Mstar-1); 
          P_exo = A_prior.rows(plag*M+pp*Mstar, plag*M+(pp+1)*Mstar-1);
          r_exo = A_exo.n_rows; c_exo = A_exo.n_cols; d_exo = A_exo.n_elem;
          
          // sample lambda
          if(pp == 1){
            prodlambda = 1.0;
          }else{
            prodlambda = as_scalar(prod(lambda2_A.submat(1,1,pp-1,1)));
          }
          dl = d_lambda + A_tau(pp,1)*d_exo;
          el = e_lambda + 0.5*A_tau(pp,1)*accu(V_end)*prodlambda;
          lambda2_A(pp,1) = R::rgamma(dl, 1/el);
          
          // sample theta
          prodlambda = as_scalar(prod(lambda2_A.submat(1,1,pp,1)));
          for(int ii=0; ii < r_exo; ii++){
            for(int jj=0; jj < c_exo; jj++){
              lambda = A_tau(pp,1) - 0.5;
              psi = A_tau(pp,1) * prodlambda;
              chi = std::pow(A_exo(ii,jj) - P_exo(ii,jj),2);
              
              res = do_rgig1(lambda, chi, psi);
              if(res<1e-7) res = 1e-7;
              if(res>1e+7) res = 1e+7;
              
              V_exo(ii,jj) = res;
            }
          }
          V_prior.rows(plag*M+pp*Mstar, plag*M+(pp+1)*Mstar-1) = V_exo;
          
          // sample tau
          if(sample_tau_bool){
            vec theta_vec = V_exo.as_col();
            
            proposal = exp(R::rnorm(0,A_tuning(pp,1)))*A_tau(pp,1);
            unif = R::runif(0,1);
            
            post_tau_prop = tau_post(proposal, prodlambda, theta_vec, 1);
            post_tau_curr = tau_post(A_tau(pp,1), prodlambda, theta_vec, 1);
            diff = post_tau_prop - post_tau_curr + std::log(proposal) - std::log(A_tau(pp,1));
            if(diff > log(unif)){
              A_tau(pp,1) = proposal;
              A_accept(pp,1) += 1;
            }
            if(irep < 0.5*burnin){
              if(A_accept(pp,1)/irep > 0.30){A_tuning(pp,1) = A_tuning(pp,1)*1.01;}
              if(A_accept(pp,1)/irep < 0.15){A_tuning(pp,1) = A_tuning(pp,1)*0.99;}
            }
          }
        }
      }
      //------------------------------------------
      // coefficients H matrix
      uvec lower_indices = trimatl_ind(size(L_draw), -1);
      vec L_vec = L_prior(lower_indices);
      
      // sample lambda
      dl = d_lambda + L_tau*v;
      el = e_lambda + 0.5*L_tau*accu(L_vec);
      lambda2_L = R::rgamma(dl, 1/el);
      
      // sample theta
      for(int ii=1; ii < M; ii++){
        for(int jj=0; jj < ii; jj++){
          lambda = L_tau - 0.5;
          psi = L_tau * lambda2_L;
          chi = std::pow(L_draw(ii,jj) - l_prior(ii,jj),2);
          
          res = do_rgig1(lambda, chi, psi);
          if(res<1e-7) res = 1e-7;
          if(res>1e+7) res = 1e+7;
          
          L_prior(ii,jj) = res;
        }
      }
      
      // sample tau
      if(sample_tau_bool){
        L_vec = L_prior(lower_indices);
        
        proposal = exp(R::rnorm(0,L_tuning))*L_tau;
        unif = R::runif(0,1);
        
        post_tau_prop = tau_post(proposal, lambda2_L, L_vec, 1);
        post_tau_curr = tau_post(L_tau, lambda2_L, L_vec, 1);
        diff= post_tau_prop - post_tau_curr + std::log(proposal) - std::log(L_tau);
        if(diff > log(unif)){
          L_tau = proposal;
          L_accept += 1;
        }
        if(irep < 0.5*burnin){
          if(L_accept/irep > 0.30){L_tuning = L_tuning*1.01;}
          if(L_accept/irep < 0.15){L_tuning = L_tuning*0.99;}
        }
      }
    }
    //-----------------------------------------------
    // Step 3: Sample covariances
    for(int mm=0; mm < M; mm++){
      vec data_sv = Em_str_draw.col(mm);
      vec cur_sv  = Sv_draw.unsafe_col(mm);  // changed to **unsafe**_col which reuses memory
      if(sv){
        const vec datastand = log(data_sv%data_sv + offset);
        mu = Sv_para(0, mm);
        phi = Sv_para(1, mm);
        sigma = Sv_para(2, mm);
        h0 = Sv_para(3, mm);
        stochvol::update_fast_sv(datastand, mu, phi, sigma, h0, cur_sv, rec, prior_spec, expert);
        Sv_para.col(mm) = arma::colvec({mu, phi, sigma});
        //Sv_draw.col(mm) = cur_sv;  // unsafe_col overwrites the original data without copying
      }else{
        a_full = a_1 + 0.5 * T;
        b_full = b_1 + 0.5 * as_scalar(data_sv.t() * data_sv);
        sig2 = 1/R::rgamma(a_full, 1/b_full);
        cur_sv.fill(sig2);
        //Sv_draw.col(mm) = log(cur_sv);
      }
    }
    //-----------------------------------------------
    // Step 4: Check Stationarity
    //gen_compMat(Cm, A_draw.rows(0,K-1), M, p);
    //cx_vec eigval; cx_mat eigvec; eig_gen(eigval, eigvec, Cm);
    //double eigval_real = abs(real(eigval)).max();
    //-----------------------------------------------
    // Step 5: STORAGE
    if(irep >= burnin){
      if((irep-burnin) % thin == 0){
        A_store.row((irep-burnin)/thin) = A_draw;
        L_store.row((irep-burnin)/thin) = L_draw;
        res_store.row((irep-burnin)/thin) = Y - X * A_draw;
        Sv_store.row((irep-burnin)/thin) = Sv_draw;
        pars_store.row((irep-burnin)/thin) = Sv_para;
        theta_store.row((irep-burnin)/thin) = V_prior;
        if(prior==1){
         shrink_store((irep-burnin)/thin,0) = shrink1;
         shrink_store((irep-burnin)/thin,1) = shrink2;
         shrink_store((irep-burnin)/thin,2) = shrink4;
        }
        if(prior==2){
          gamma_store.row((irep-burnin)/thin) = gamma;
          omega_store.row((irep-burnin)/thin) = omega;
        }
        if(prior==3){
          lambda2_Lmat(0,0) = lambda2_L;
          lambda2_store.row((irep-burnin)/thin) = join_rows(lambda2_A,lambda2_Lmat);
          L_taumat(0,0) = L_tau;
          tau_store.row((irep-burnin)/thin) = join_rows(A_tau,L_taumat);
        }
      }
    }
    // check user interruption
    if(irep % 200 == 0)
      Rcpp::checkUserInterrupt();
  } // END MCMC LOOP
  //---------------------------------------------------------------------------------------------
  return List::create(Named("Y") = Y, 
                      Named("X") = X, 
                      Named("A_store") = A_store, 
                      Named("L_store") = L_store, 
                      Named("Sv_store") = Sv_store, 
                      Named("shrink_store") = shrink_store,
                      Named("gamma_store") = gamma_store, 
                      Named("omega_store") = omega_store, 
                      Named("theta_store") = theta_store, 
                      Named("lambda2_store") = lambda2_store, 
                      Named("tau_store") = tau_store, 
                      Named("pars_store") = pars_store, 
                      Named("res_store") = res_store);
}
