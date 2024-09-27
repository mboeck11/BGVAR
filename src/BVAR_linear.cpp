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
                 arma::uvec lags,
                 int draws, 
                 int burnin,
                 int thin,
                 bool cons, 
                 bool trend,
                 bool sv,
                 int prior,
                 Rcpp::List hyperparam,
                 Rcpp::List setting_store) {
  //----------------------------------------------------------------------------------------------------------------------
  // CONSTRUCT DATA
  //----------------------------------------------------------------------------------------------------------------------
  int plag = lags(0);
  int plagstar = lags(1);
  int pmax  = max(lags);
  
  int Traw = Yraw.n_rows;
  int M = Yraw.n_cols;
  int K = M*plag;
  int Mstar = Wraw.n_cols;
  int Kstar = Mstar*(plagstar+1);
  
  bool texo = false; int Mex=0;
  if(Exraw.n_elem != 1){
    texo = true; Mex = Exraw.n_cols; 
  }
  
  mat Xraw = mlag(Yraw,plag,Traw,M);
  mat X0   = Xraw.submat(pmax,0,Traw-1,K-1);
  mat Y    = Yraw.submat(pmax,0,Traw-1,M-1);
  double T = X0.n_rows;
  mat Wall = join_rows(Wraw,mlag(Wraw,plagstar,Traw,Mstar));
  mat W0   = Wall.submat(pmax,0,Traw-1,Kstar-1);
  mat X    = join_rows(X0,W0);
  if(texo){
    mat E0 = Exraw.submat(pmax,0,Traw-1,Mex-1);
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
  int n = K*M;
  int nstar = Kstar*M;
  
  //----------------------------------------------------------------------------------------------------------------------
  // HYPERPARAMETERS
  //----------------------------------------------------------------------------------------------------------------------
  // prior == 1: SIMS
  double lambda1 = hyperparam["lambda1"];
  double lambda2 = hyperparam["lambda2"];
  double lambda3 = hyperparam["lambda3"];
  double lambda4 = hyperparam["lambda4"];
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
  // STORE SETTINGS
  //----------------------------------------------------------------------------------------------------------------------
  const bool save_shrink_MN   = setting_store["shrink_MN"];
  const bool save_shrink_SSVS = setting_store["shrink_SSVS"];
  const bool save_shrink_NG   = setting_store["shrink_NG"];
  const bool save_shrink_HS   = setting_store["shrink_HS"];
  const bool save_vola_pars   = setting_store["vola_pars"];
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
  mat L_drawinv(M,M); L_drawinv.eye();
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
  double accept1 = 0, accept2 = 0, accept3 = 0;
  double scale1 = 0.43, scale2 = 0.43, scale3 = 0.43;
  vec sigmas(M+Mstar, fill::zeros);
  for(int i=0; i < M; i++){
    mat Y_ = Yraw.col(i);
    sigmas(i) = get_ar(Y_,plag);
  }
  for(int i=0; i < Mstar; i++){
    mat W_ = Wraw.col(i);
    sigmas(M+i) = get_ar(W_,plagstar);
  }
  if(prior==1) get_Vminnesota(V_prior, sigmas, lambda1, lambda2, lambda3, lambda4, cons, Mstar, plag, plagstar, trend);
  
  // initialize stuff for MN prior
  mat V_prop1(k,M), V_prop2(k,M), V_prop3(k,M); 
  double lambda_prop1, lambda_prop2, lambda_prop3;
  double post_old1=0.0, post_prop1=0.0, post_old2=0.0, post_prop2=0.0, post_old3=0.0, post_prop3=0.0;
  for(int i=0; i<k; i++){
    for(int j=0; j<M; j++){
      post_old1 = post_old1 + R::dnorm(A_draw(i,j),A_prior(i,j),  std::sqrt(V_prior(i,j)),true);
      post_old2 = post_old2 + R::dnorm(A_draw(i,j),A_prior(i,j),  std::sqrt(V_prior(i,j)),true);
      post_old3 = post_old3 + R::dnorm(A_draw(i,j), A_prior(i,j), std::sqrt(V_prior(i,j)), true);
    }
  }
  post_old1 = post_old1 + R::dgamma(lambda1,0.01,1/0.01,true) + log(lambda1); // add prior - shape scale parameterization!!!! + correction term
  post_old2 = post_old2 + R::dgamma(lambda2,0.01,1/0.01,true) + log(lambda2); // add prior - shape scale parameterization!!!! + correction term
  post_old3 = post_old3 + R::dgamma(lambda3,0.01,1/0.01,true) + log(lambda4); // add prior - shape scale parameterization!!!! + correction term
  
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
  mat lambda2_A(pmax+1,2,fill::zeros);
  mat A_tau(pmax+1,2); A_tau.fill(tau_theta); A_tau(0,0)=0;
  mat A_tuning(pmax+1,2); A_tuning.fill(0.43);
  mat A_accept(pmax+1,2, fill::zeros);
  // initialize stuff for NG prior
  mat A_con, V_con, P_con, A_end, V_end, P_end, A_exo, V_exo, P_exo;
  double prodlambda, dl, el, lambda, chi, psi, res;
  double unif, proposal, post_tau_prop, post_tau_curr, diff;
  
  // HS stuff
  vec lambda_A_endo(n, fill::ones), nu_A_endo(n, fill::ones), 
      lambda_A_exo(nstar, fill::ones), nu_A_exo(nstar, fill::ones),
      lambda_L(v, fill::ones), nu_L(v, fill::ones);
  double tau_A_endo = 1, zeta_A_endo = 1, tau_A_exo = 1, zeta_A_exo = 1, tau_L = 1, zeta_L = 1;
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
  mat lambda2_Lmat(pmax+1,1, fill::zeros);
  mat L_taumat(pmax+1,1, fill::zeros); 
  mat L_accmat(pmax+1,1, fill::zeros);
  mat L_tunmat(pmax+1,1, fill::zeros);
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
  const double offset = 1e-40;  // maybe want to change to 1e-40 or so to be on the safe side? I have got random NAs because of log(0) for real data inputs
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
  arma::cube A_store(k,M,thindraws);
  arma::cube L_store(M,M,thindraws);
  arma::cube res_store(T,M,thindraws);
  // SV
  arma::cube Sv_store(T,M,thindraws);
  arma::ivec size_of_cube1 = {0, 0, 0};
  arma::ivec size_of_cube2 = {0, 0, 0};
  if(save_vola_pars == true){
    size_of_cube1 = {4, M, thindraws};
  }
  arma::cube pars_store(size_of_cube1(0), size_of_cube1(1), size_of_cube1(2));
  // MN
  if(save_shrink_MN == true){
    size_of_cube2 = {3, 1, thindraws};
  }
  arma::cube lambda_store(size_of_cube2(0), size_of_cube2(1), size_of_cube2(2));
  // SSVS
  size_of_cube1 = {0, 0, 0};
  size_of_cube2 = {0, 0, 0};
  if(save_shrink_SSVS == true){
    size_of_cube1 = {k, M, thindraws};
    size_of_cube2 = {M, M, thindraws};
  }
  arma::cube gamma_store(size_of_cube1(0), size_of_cube1(1), size_of_cube1(2));
  arma::cube omega_store(size_of_cube2(0), size_of_cube2(1), size_of_cube2(2));
  // NG
  if(save_shrink_NG){
    size_of_cube1 = {k, M, thindraws};
    size_of_cube2 = {pmax+1, 3, thindraws};
  }
  arma::cube theta_store(size_of_cube1(0), size_of_cube1(1), size_of_cube1(2));
  arma::cube lambda2_store(size_of_cube2(0), size_of_cube2(1), size_of_cube2(2));
  arma::cube tau_store(size_of_cube2(0), size_of_cube2(1), size_of_cube2(2));
  // HS
  arma::ivec size_of_mat1 = {0, 0}, size_of_mat2 = {0, 0}, size_of_mat3 = {0, 0}, size_of_mat4 = {0, 0};
  if(save_shrink_HS){
    size_of_mat1 = {n, thindraws};
    size_of_mat2 = {nstar, thindraws};
    size_of_mat3 = {v, thindraws};
    size_of_mat4 = {1, thindraws};
  }
  arma::mat lambda_A_endo_store(size_of_mat1(0), size_of_mat1(1));
  arma::mat lambda_A_exo_store(size_of_mat2(0), size_of_mat2(1));
  arma::mat lambda_L_store(size_of_mat3(0), size_of_mat3(1));
  arma::mat nu_A_endo_store(size_of_mat1(0), size_of_mat1(1));
  arma::mat nu_A_exo_store(size_of_mat2(0), size_of_mat2(1));
  arma::mat nu_L_store(size_of_mat3(0), size_of_mat3(1));
  arma::mat tau_A_endo_store(size_of_mat4(0), size_of_mat4(1));
  arma::mat tau_A_exo_store(size_of_mat4(0), size_of_mat4(1));
  arma::mat tau_L_store(size_of_mat4(0), size_of_mat4(1));
  arma::mat zeta_A_endo_store(size_of_mat4(0), size_of_mat4(1));
  arma::mat zeta_A_exo_store(size_of_mat4(0), size_of_mat4(1));
  arma::mat zeta_L_store(size_of_mat4(0), size_of_mat4(1));
  //---------------------------------------------------------------------------------------------
  // MCMC LOOP
  //---------------------------------------------------------------------------------------------
  for(int irep = 0; irep < ntot; irep++){
    // Step 1: Sample coefficients
    // Step 1a: Sample coefficients in A
    for(int mm = 0; mm < M; mm++){ // estimate equation-by-equation
      mat A_0    = A_draw; A_0.col(mm).fill(0);
      mat Linv_0 = L_drawinv.rows(mm,M-1);
      mat S_0    = exp(-0.5*Sv_draw.cols(mm,M-1));
      mat zmat   = (Y - X * A_0) * Linv_0.t() ;
      vec ztilde = vectorise(zmat) % vectorise(S_0);
      mat xtilde = kron(Linv_0.col(mm), X) % repmat(vectorise(S_0),1,k);
      mat Vinv_m = diagmat(1/V_prior.col(mm));
      colvec a_m = A_prior.col(mm);
      
      mat V_p = (xtilde.t() * xtilde + Vinv_m).i();
      mat A_p = V_p * (xtilde.t() * ztilde + Vinv_m * a_m);
      
      colvec rand_normal(k);
      for(int i=0; i<k; i++){
        rand_normal(i) = R::rnorm(0,1);
      }
      mat V_p_chol_lower = robust_chol(V_p);
      /*
      bool chol_success = chol(V_p_chol_lower, V_p, "lower");
      // Fall back on Rs chol if armadillo fails (it suppports pivoting)
      if(chol_success == false){
        NumericMatrix tmp = Rchol(V_p, true);
        int d = V_p.n_cols;
        mat cholV_tmp = mat(tmp.begin(), d, d, false);
        uvec piv = sort_index(as<vec>(tmp.attr("pivot")));
        V_p_chol_lower = cholV_tmp.cols(piv);
        V_p_chol_lower = V_p_chol_lower.t();
      }
       */
      colvec A_m = A_p + V_p_chol_lower*rand_normal;
      
      A_draw.col(mm) = A_m;
      Em_draw.col(mm) = Y.col(mm) - X * A_m;
    }
    // Step 1b: Sample coefficients in L
    for(int mm = 1; mm < M; mm++){ 
      mat eps_m = Em_draw.col(mm);
      mat S_m   = exp(-0.5*Sv_draw.col(mm));
      mat eps_x = Em_draw.cols(0,mm-1);
      eps_m = eps_m % S_m;
      eps_x = eps_x % repmat(S_m,1,mm);
      mat Vinv_m = diagmat(1/L_prior.submat(mm,0,mm,mm-1));
      colvec a_m = l_prior.submat(mm,0,mm,mm-1).t();
      
      mat V_p = (eps_x.t() * eps_x + Vinv_m).i();
      mat A_p = V_p * (eps_x.t() * eps_m + Vinv_m * a_m);
      
      colvec rand_normal(mm);
      for(int i=0; i< mm; i++){
        rand_normal(i) = R::rnorm(0,1);
      }
      mat V_p_chol_lower = robust_chol(V_p);
      /*
      bool chol_success = chol(V_p_chol_lower, V_p,"lower");
      // Fall back on Rs chol if armadillo fails (it suppports pivoting)
      if(chol_success == false){
        NumericMatrix tmp = Rchol(V_p, true);
        int d = V_p.n_cols;
        mat cholV_tmp = mat(tmp.begin(), d, d, false);
        uvec piv = sort_index(as<vec>(tmp.attr("pivot")));
        V_p_chol_lower = cholV_tmp.cols(piv);
        V_p_chol_lower = V_p_chol_lower.t();
      }
       */
      colvec L_m = A_p + V_p_chol_lower*rand_normal;
      
      L_draw.submat(mm,0,mm,mm-1) = L_m;
    }
    L_drawinv = L_draw.i();
    Em_str_draw = Y * L_drawinv.t() - X * A_draw * L_drawinv.t();
    /*
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
          NumericMatrix tmp = Rchol(V_p, true);
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
          NumericMatrix tmp = Rchol(V_p, true);
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
     */
    //-----------------------------------------------
    // Step 2: different prior setups
    // SIMS
    if(prior == 1){
      // first shrinkage parameter (own lags)
      lambda_prop1 = exp(R::rnorm(0,scale1))*lambda1;
      get_Vminnesota(V_prop1, sigmas, lambda_prop1, lambda2, lambda3, lambda4, cons, Mstar, plag, plagstar, trend);
      // likelihood of each coefficient
      for(int i=0; i<k; i++){
        for(int j=0; j<M; j++){
          post_prop1 = post_prop1 + R::dnorm(A_draw(i,j), A_prior(i,j), std::sqrt(V_prop1(i,j)), true);
        }
      }
      // total likelihood 
      post_prop1 = post_prop1 + R::dgamma(lambda_prop1,0.01,1/0.01,true) + log(lambda_prop1);  // add prior - shape scale parameterization!!!! + correction term
      if((post_prop1-post_old1) > log(R::runif(0,1))){
        lambda1 = lambda_prop1;
        V_prior = V_prop1;
        post_old1 = post_prop1;
        accept1 += 1;
      }
      
      // second shrinkage parameter (cross equations)
      lambda_prop2 = exp(R::rnorm(0,scale2))*lambda2;
      get_Vminnesota(V_prop2, sigmas, lambda1, lambda_prop2, lambda3, lambda4, cons, Mstar, plag, plagstar, trend);
      // likelihood of each coefficient
      for(int i=0; i<k; i++){
        for(int j=0; j<M; j++){
          post_prop2 = post_prop2 + R::dnorm(A_draw(i,j), A_prior(i,j), std::sqrt(V_prop2(i,j)), true);
        }
      }
      // total likelihood 
      post_prop2 = post_prop2 + R::dgamma(lambda_prop2,0.01,1/0.01,true) + log(lambda_prop2);  // add prior - shape scale parameterization!!!! + correction term
      if((post_prop2-post_old2) > log(R::runif(0,1))){
        lambda2 = lambda_prop2;
        V_prior = V_prop2;
        post_old2 = post_prop2; 
        accept2 += 1;
      }
      
      // third shrinkage parameter (weakly exogenous)
      lambda_prop3 = exp(R::rnorm(0,scale3))*lambda3;
      get_Vminnesota(V_prop3, sigmas, lambda1, lambda2, lambda3, lambda_prop3, cons, Mstar, plag, plagstar, trend);
      // likelihood of each coefficient
      for(int i=0; i<k; i++){
        for(int j=0; j<M; j++){
          post_prop3 = post_prop3 + R::dnorm(A_draw(i,j), A_prior(i,j), std::sqrt(V_prop3(i,j)), true);
        }
      }
      // total likelihood 
      post_prop3 = post_prop3 + R::dgamma(lambda_prop3,0.01,1/0.01,true) + log(lambda_prop3);  // add prior - shape scale parameterization!!!! + correction term
      if((post_prop3-post_old3) > log(R::runif(0,1))){
        lambda3 = lambda_prop3;
        V_prior = V_prop3;
        post_old3 = post_prop3;
        accept3 += 1;
      }
      
      if((irep+1) < 0.5*burnin){
        if(accept1/(irep+1) > 0.30){scale1 *= 1.01;}
        if(accept1/(irep+1) < 0.15){scale1 *= 0.99;}
        if(accept2/(irep+1) > 0.30){scale2 *= 1.01;}
        if(accept2/(irep+1) < 0.15){scale2 *= 0.99;}
        if(accept3/(irep+1) > 0.30){scale3 *= 1.01;}
        if(accept3/(irep+1) < 0.15){scale3 *= 0.99;}
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
      // coefficients A matrix - weakly exogenous
      for(int pp=0; pp < (plagstar+1); pp++){
        //-------------------------------------------------------------------
        // weakly exogenous
        A_exo = A_draw.rows(plag*M+pp*Mstar, plag*M+(pp+1)*Mstar-1); 
        V_exo = V_prior.rows(plag*M+pp*Mstar, plag*M+(pp+1)*Mstar-1); 
        P_exo = A_prior.rows(plag*M+pp*Mstar, plag*M+(pp+1)*Mstar-1);
        
        // sample lambda
        if(pp == 0){
          prodlambda = 1.0;
        }else{
          prodlambda = as_scalar(prod(lambda2_A.submat(0,1,pp-1,1)));
        }
        dl = d_lambda + A_tau(pp,1)*M*Mstar;
        el = e_lambda + 0.5*A_tau(pp,1)*arma::accu(V_exo)*prodlambda;
        lambda2_A(pp,1) = R::rgamma(dl, 1/el);
        
        // sample theta
        prodlambda = as_scalar(prod(lambda2_A.submat(0,1,pp,1)));
        for(int ii=0; ii < Mstar; ii++){
          for(int mm=0; mm < M; mm++){
            lambda = A_tau(pp,1) - 0.5;
            psi = A_tau(pp,1) * prodlambda;
            chi = std::pow(A_exo(ii,mm) - P_exo(ii,mm),2);
            
            res = do_rgig1(lambda, chi, psi);
            if(res<1e-7) res = 1e-7;
            if(res>1e+7) res = 1e+7;
            
            V_exo(ii,mm) = res;
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
      // coefficients A matrix
      for(int pp=0; pp < plag; pp++){
        A_end = A_draw.rows(pp*M, (pp+1)*M-1); 
        V_end = V_prior.rows(pp*M, (pp+1)*M-1); 
        P_end = A_prior.rows(pp*M, (pp+1)*M-1);
        
        // sample lambda
        if(pp == 0){
          prodlambda = 1.0;
        }else{
          prodlambda = as_scalar(prod(lambda2_A.submat(1,0,pp,0)));
        }
        dl = d_lambda + A_tau(pp+1,0)*std::pow(M,2);
        el = e_lambda + 0.5*A_tau(pp+1,0)*arma::accu(V_end)*prodlambda;
        lambda2_A(pp+1,0) = R::rgamma(dl, 1/el);
        
        // sample theta
        prodlambda = as_scalar(prod(lambda2_A.submat(1,0,pp+1,0)));
        for(int ii=0; ii < M; ii++){
          for(int mm=0; mm < M; mm++){
            lambda = A_tau(pp+1,0) - 0.5;
            psi = A_tau(pp+1,0) * prodlambda;
            chi = std::pow(A_end(ii,mm)-P_end(ii,mm),2);
            
            res = do_rgig1(lambda, chi, psi);
            if(res<1e-7) res = 1e-7;
            if(res>1e+7) res = 1e+7;
            
            V_end(ii,mm) = res;
          }
        }
        V_prior.rows(pp*M, (pp+1)*M-1) = V_end;
        
        // sample tau
        if(sample_tau_bool){
          vec theta_vec = V_end.as_col();
          
          proposal = exp(R::rnorm(0,A_tuning(pp+1,0)))*A_tau(pp+1,0);
          unif = R::runif(0,1);
          
          post_tau_prop = tau_post(proposal, prodlambda, theta_vec, 1);
          post_tau_curr = tau_post(A_tau(pp+1,0), prodlambda, theta_vec, 1);
          diff = post_tau_prop - post_tau_curr + std::log(proposal) - std::log(A_tau(pp+1,0));
          if(diff > log(unif)){
            A_tau(pp+1,0) = proposal;
            A_accept(pp+1,0) += 1;
          }
          if(irep < 0.5*burnin){
            if(A_accept(pp+1,0)/irep > 0.30){A_tuning(pp+1,0) = A_tuning(pp+1,0)*1.01;}
            if(A_accept(pp+1,0)/irep < 0.15){A_tuning(pp+1,0) = A_tuning(pp+1,0)*0.99;}
          }
        }
      }
      //------------------------------------------
      // coefficients H matrix
      uvec lower_indices = trimatl_ind(size(L_draw), -1);
      vec L_vec = L_prior(lower_indices);
      
      // sample lambda
      dl = d_lambda + L_tau*v;
      el = e_lambda + 0.5*L_tau*arma::accu(L_vec);
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
    // HS
    if(prior == 4){
      //------------------------------------------
      // coefficients H matrix
      uvec lower_indices = trimatl_ind(size(L_draw), -1);
      // sample local shrinkage parameter
      for(int vv=0; vv < v; vv++){
        lambda_L(vv) = 1/R::rgamma(1, 1/(1 / nu_L(vv) + 0.5*pow(L_draw(lower_indices(vv)),2) / tau_L));
        nu_L(vv)     = 1/R::rgamma(1, 1/(1 + 1/lambda_L(vv)));
      }
      // sample global shrinkage parameter
      tau_L  = 1/R::rgamma((v+1)/2, 1/(1/zeta_L + 0.5*sum(pow(L_draw(lower_indices),2)/lambda_L)));
      zeta_L = 1/R::rgamma(1, 1/(1 + 1 / tau_L));
      // update prior VCV
      L_prior(lower_indices) = tau_L * lambda_L;
      
      //------------------------------------------
      // coefficients A matrix - endogenous
      A_end = A_draw.rows(0, plag*M-1);
      // sample local shrinkage parameter
      for(int nn=0; nn < n; nn++){
        lambda_A_endo(nn) = 1.0 / R::rgamma(1, 1/(1 / nu_A_endo(nn) + 0.5*pow(A_end(nn),2) / tau_A_endo));
        nu_A_endo(nn)     = 1.0 / R::rgamma(1, 1/(1 + 1/lambda_A_endo(nn)));
      }
      // sample global shrinkage parameter
      tau_A_endo  = 1.0/R::rgamma((n+1)/2, 1/(1/zeta_A_endo + 0.5*sum(pow(vectorise(A_end),2)/lambda_A_endo)));
      zeta_A_endo = 1.0/R::rgamma(1, 1 + 1/(1 / tau_A_endo));
      // update prior VCV
      V_prior.rows(0, plag*M-1) = reshape(tau_A_endo * lambda_A_endo, plag*M, M);
      
      //------------------------------------------
      // coefficients A matrix - exogenous
      A_exo = A_draw.rows(plag*M, plag*M+(plagstar+1)*Mstar-1);
      // sample local shrinkage parameter
      for(int nn=0; nn < nstar; nn++){
        lambda_A_exo(nn) = 1/R::rgamma(1, 1/(1 / nu_A_exo(nn) + 0.5*pow(A_exo(nn),2) / tau_A_exo));
        nu_A_exo(nn)     = 1/R::rgamma(1, 1/(1 + 1/lambda_A_exo(nn)));
      }
      // sample global shrinkage parameter
      tau_A_exo  = 1/R::rgamma((nstar+1)/2, 1/(1/zeta_A_exo + 0.5*sum(pow(vectorise(A_exo),2)/lambda_A_exo)));
      zeta_A_exo = 1/R::rgamma(1, 1/(1 + 1 / tau_A_exo));
      // update prior VCV
      V_prior.rows(plag*M, plag*M+(plagstar+1)*Mstar-1) = reshape(tau_A_exo * lambda_A_exo, (plagstar+1)*Mstar, M);
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
        Sv_para.col(mm) = arma::colvec({mu, phi, sigma, h0});
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
        A_store.slice((irep-burnin)/thin) = A_draw;
        L_store.slice((irep-burnin)/thin) = L_draw;
        res_store.slice((irep-burnin)/thin) = Y - X * A_draw;
        Sv_store.slice((irep-burnin)/thin) = Sv_draw;
        if(save_vola_pars == true){
          pars_store.slice((irep-burnin)/thin) = Sv_para;
        }
        if(save_shrink_MN == true){
         lambda_store(0,0,(irep-burnin)/thin) = lambda1;
         lambda_store(1,0,(irep-burnin)/thin) = lambda2;
         lambda_store(2,0,(irep-burnin)/thin) = lambda3;
        }
        if(save_shrink_SSVS == true){
          gamma_store.slice((irep-burnin)/thin) = gamma;
          omega_store.slice((irep-burnin)/thin) = omega;
        }
        if(save_shrink_NG == true){
          theta_store.slice((irep-burnin)/thin) = V_prior;
          lambda2_Lmat(0,0) = lambda2_L;
          lambda2_store.slice((irep-burnin)/thin) = join_rows(lambda2_A,lambda2_Lmat);
          L_taumat(0,0) = L_tau;
          tau_store.slice((irep-burnin)/thin) = join_rows(A_tau,L_taumat);
        }
        if(save_shrink_HS == true){
          lambda_A_endo_store.col((irep-burnin)/thin) = lambda_A_endo;
          lambda_A_exo_store.col((irep-burnin)/thin)  = lambda_A_exo;
          lambda_L_store.col((irep-burnin)/thin)      = lambda_L;
          nu_A_endo_store.col((irep-burnin)/thin)     = nu_A_endo;
          nu_A_exo_store.col((irep-burnin)/thin)      = nu_A_exo;
          nu_L_store.col((irep-burnin)/thin)          = nu_L;
          tau_A_endo_store.col((irep-burnin)/thin)    = tau_A_endo;
          tau_A_exo_store.col((irep-burnin)/thin)     = tau_A_exo;
          tau_L_store.col((irep-burnin)/thin)         = tau_L;
          zeta_A_endo_store.col((irep-burnin)/thin)   = zeta_A_endo;
          zeta_A_exo_store.col((irep-burnin)/thin)    = zeta_A_exo;
          zeta_L_store.col((irep-burnin)/thin)        = zeta_L;
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
                      Named("res_store") = res_store,
                      Named("pars_store") = pars_store, 
                      Named("MN") = List::create(
                        Named("lambda_store") = lambda_store
                      ),
                      Named("SSVS") = List::create(
                        Named("gamma_store") = gamma_store, 
                        Named("omega_store") = omega_store
                      ),
                      Named("NG") = List::create(
                        Named("theta_store") = theta_store, 
                        Named("lambda2_store") = lambda2_store, 
                        Named("tau_store") = tau_store
                      ),
                      Named("HS") = List::create(
                        Named("lambda_A_endo_store") = lambda_A_endo_store,
                        Named("lambda_A_exo_store") = lambda_A_exo_store,
                        Named("lambda_L_store") = lambda_L_store,
                        Named("nu_A_endo_store") = nu_A_endo_store,
                        Named("nu_A_exo_store") = nu_A_exo_store,
                        Named("nu_L_store") = nu_L_store,
                        Named("tau_A_exo_store") = tau_A_exo_store,
                        Named("tau_A_endo_store") = tau_A_endo_store,
                        Named("tau_L_store") = tau_L_store,
                        Named("zeta_A_endo_store") = zeta_A_endo_store,
                        Named("zeta_A_exo_store") = zeta_A_exo_store,
                        Named("zeta_L_store") = zeta_L_store)
                      );
}
