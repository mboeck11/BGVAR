#include <RcppArmadillo.h>
#include <stochvol.h>
#include <stdlib.h>
#include "helper.h"
#include "sample_parameters.h"
#include "do_rgig1.h"

using namespace Rcpp;
using namespace arma;

//' @name BVAR_linear
//' @noRd
//[[Rcpp::interfaces(r, cpp)]]
//[[Rcpp::export]]
List BVAR_linear(const SEXP Y_in, const SEXP W_in, const SEXP p_in,
                 const SEXP draws_in, const SEXP burnin_in,
                 const SEXP cons_in, const SEXP trend_in, const SEXP sv_in, const SEXP thin_in,
                 const SEXP prior_in, const SEXP hyperparam_in, const SEXP Ex_in) {
  //----------------------------------------------------------------------------------------------------------------------
  // CONSTRUCT DATA
  //----------------------------------------------------------------------------------------------------------------------
  NumericMatrix Yr(Y_in);
  int Traw = Yr.nrow(), M = Yr.ncol();
  mat Yraw(Yr.begin(), Traw, M, false);
  const int p = as<int>(p_in);
  int K = M*p;

  bool exo = false;
  mat Wraw; NumericMatrix Wr; int Mstar=0; int Kstar=0;
  if(W_in != R_NilValue) {
    exo = true; NumericMatrix Wr(W_in); Mstar = Wr.ncol(); Kstar = Mstar*(p+1);
    Wraw = mat(Wr.begin(), Traw, Mstar, false);
  }
  
  bool texo = false;
  mat Exraw; NumericMatrix Er; int Mex=0;
  if(Ex_in != R_NilValue){
    texo = true; NumericMatrix Er(Ex_in); Mex = Er.ncol(); 
    Exraw = mat(Er.begin(), Traw, Mex, false);
  }
  
  mat Xraw = mlag(Yraw,p,Traw,M);
  mat X0 = Xraw.submat(p,0,Traw-1,K-1);
  mat X = X0;
  mat Y = Yraw.submat(p,0,Traw-1,M-1);
  double T = X0.n_rows;
  if(exo){
    mat Wall = join_rows(Wraw,mlag(Wraw,p,Traw,Mstar));
    mat W0 = Wall.submat(p,0,Traw-1,Kstar-1);
    X = join_rows(X0,W0);
  }
  if(texo){
    mat E0 = Exraw.submat(p,0,Traw-1,Mex-1);
    X = join_rows(X,E0);
  }
  const bool cons = as<bool>(cons_in);
  if(cons){
    X = join_rows(X,colvec(T,fill::ones));
  }
  const bool trend = as<bool>(trend_in);
  if(trend){
    vec trendvec(T); for(int tt=0; tt<T; tt++) {trendvec(tt)=tt+1;}
    X = join_rows(X,trendvec);
  }

  int k = X.n_cols;
  int v = (M*(M-1))/2;
  //----------------------------------------------------------------------------------------------------------------------
  // HYPERPARAMETERS
  //----------------------------------------------------------------------------------------------------------------------
  List hyperparam(clone(hyperparam_in));
  const int prior = as<int>(prior_in);
  const bool sv = as<int>(sv_in);
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
  const double a_start = hyperparam["a_start"];
  const bool sample_A = hyperparam["sample_A"];
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
  List coeflist;
  mat A_draw = A_OLS;
  cube SIGMA(M,M,T);
  for(int tt=0; tt < T; tt++){
    SIGMA.slice(tt) = SIGMA_OLS;
  }
  mat Em_draw = E_OLS; mat Em_str_draw = E_OLS;
  mat L_draw(M,M); L_draw.eye();
  mat Cm(K,K, fill::zeros); gen_compMat(Cm, A_OLS.rows(0,K-1), M, p);
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
    sigmas(i) = get_ar(Y_,p);
  }
  if(exo){
    for(int i=0; i < Mstar; i++){
      mat W_ = Wraw.col(i);
      sigmas(M+i) = get_ar(W_,p);
    }
  }
  if(prior==1) get_Vminnesota(V_prior, sigmas, shrink1, shrink2, shrink3, shrink4, cons, Mstar, p, trend);
  
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
  mat lambda2_A(p+1,2,fill::zeros);
  mat A_tau(p+1,2); A_tau.fill(a_start);
  mat A_tuning(p+1,2); A_tuning.fill(0.43);
  mat A_accept(p+1,2, fill::zeros);
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
  double L_tau = a_start;
  double L_tuning = 0.43;
  double L_accept = 0;
  mat lambda2_Lmat(p+1,1, fill::zeros);
  mat L_taumat(p+1,1, fill::zeros); 
  mat L_accmat(p+1,1, fill::zeros);
  mat L_tunmat(p+1,1, fill::zeros);
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
  double h0 = -10;
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
  //---------------------------------------------------------------------------------------------------------------
  // SAMPLER MISCELLANEOUS
  //---------------------------------------------------------------------------------------------------------------
  const int burnin = as<int>(burnin_in);
  const int draws  = as<int>(draws_in);
  const int ntot   = burnin + draws;

  // thinning parameters
  const int thin = as<int>(thin_in);
  const int thindraws = draws/thin;
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
  cube lambda2_store(thindraws,p+1,3, fill::zeros);
  cube tau_store(thindraws,p+1,3, fill::zeros);
  //cube accept_store(thindraws,p+1,3, fill::zeros);
  //cube tuning_store(thindraws,p+1,3, fill::zeros);
  //---------------------------------------------------------------------------------------------
  // MCMC LOOP
  //---------------------------------------------------------------------------------------------
  for(int irep = 0; irep < ntot; irep++){
    // Step 1: Sample coefficients
    sample_arcoefs(A_draw, L_draw, Em_draw, Em_str_draw, Y, X, Sv_draw, A_prior, V_prior, l_prior, L_prior);
    //-----------------------------------------------
    // Step 2: different prior setups
    // SIMS
    if(prior == 1){
      mat postval_prop1(k,M, fill::zeros); mat postval_old1(k,M, fill::zeros);
      mat postval_prop2(k,M, fill::zeros); mat postval_old2(k,M, fill::zeros);
      mat postval_prop4(k,M, fill::zeros); mat postval_old4(k,M, fill::zeros);
      // first shrinkage parameter (own lags)
      double shrink_prop1 = exp(R::rnorm(0,scale1))*shrink1;
      mat V_prop1(k,M); get_Vminnesota(V_prop1, sigmas, shrink_prop1, shrink2, shrink3, shrink4, cons, Mstar, p, trend);
      //get_shrink(V_prior, V_prop1, shrink1, shrink_prop1, A_draw, A_prior, accept1, scale1, irep, burnin);
      
      // likelihood of each coefficient
      for(int i=0; i<k; i++){
        for(int j=0; j<M; j++){
          postval_prop1(i,j) = R::dnorm(A_draw(i,j),A_prior(i,j),std::sqrt(V_prop1(i,j)),true);
          postval_old1(i,j) = R::dnorm(A_draw(i,j),A_prior(i,j),std::sqrt(V_prior(i,j)),true);
        }
      }
      // total likelihood 
      double post_prop1 = accu(postval_prop1) + R::dgamma(shrink_prop1,0.01,1/0.01,true) + log(shrink_prop1);  // add prior - shape scale parameterization!!!! + correction term
      double post_old1 = accu(postval_old1) + R::dgamma(shrink1,0.01,1/0.01,true) + log(shrink1); // add prior - shape scale parameterization!!!! + correction term
      if((post_prop1-post_old1) > log(R::runif(0,1))){
        shrink1 = shrink_prop1;
        V_prior = V_prop1;
        accept1 += 1;
      }
      
      // second shrinkage parameter (cross equations)
      double shrink_prop2 = exp(R::rnorm(0,scale2))*shrink2;
      mat V_prop2(k,M); get_Vminnesota(V_prop2, sigmas, shrink1, shrink_prop2, shrink3, shrink4, cons, Mstar, p, trend);
      //get_shrink(V_prior, V_prop2, shrink2, shrink_prop2, A_draw, A_prior, accept2, scale2, irep, burnin);
      
      // likelihood of each coefficient
      for(int i=0; i<k; i++){
        for(int j=0; j<M; j++){
          postval_prop2(i,j) = R::dnorm(A_draw(i,j),A_prior(i,j),std::sqrt(V_prop2(i,j)),true);
          postval_old2(i,j) = R::dnorm(A_draw(i,j),A_prior(i,j),std::sqrt(V_prior(i,j)),true);
        }
      }
      // total likelihood 
      double post_prop2 = accu(postval_prop2) + R::dgamma(shrink_prop2,0.01,1/0.01,true) + log(shrink_prop2);  // add prior - shape scale parameterization!!!! + correction term
      double post_old2 = accu(postval_old2) + R::dgamma(shrink2,0.01,1/0.01,true) + log(shrink2); // add prior - shape scale parameterization!!!! + correction term
      if((post_prop2-post_old2) > log(R::runif(0,1))){
        shrink2 = shrink_prop2;
        V_prior = V_prop2;
        accept2 += 1;
      }
      
      // fourth shrinkage parameter (weakly exogenous)
      double shrink_prop4 = exp(R::rnorm(0,scale4))*shrink4;
      mat V_prop4(k,M); get_Vminnesota(V_prop4, sigmas, shrink1, shrink2, shrink3, shrink_prop4, cons, Mstar, p, trend);
      //get_shrink(V_prior, V_prop4, shrink4, shrink_prop4, A_draw, A_prior, accept4, scale4, irep, burnin);
      
      // likelihood of each coefficient
      for(int i=0; i<k; i++){
        for(int j=0; j<M; j++){
          postval_prop4(i,j) = R::dnorm(A_draw(i,j),A_prior(i,j),std::sqrt(V_prop4(i,j)),true);
          postval_old4(i,j) = R::dnorm(A_draw(i,j),A_prior(i,j),std::sqrt(V_prior(i,j)),true);
        }
      }
      // total likelihood 
      double post_prop4 = accu(postval_prop4) + R::dgamma(shrink_prop4,0.01,1/0.01,true) + log(shrink_prop4);  // add prior - shape scale parameterization!!!! + correction term
      double post_old4 = accu(postval_old4) + R::dgamma(shrink4,0.01,1/0.01,true) + log(shrink4); // add prior - shape scale parameterization!!!! + correction term
      if((post_prop4-post_old4) > log(R::runif(0,1))){
        shrink4 = shrink_prop4;
        V_prior = V_prop4;
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
          double u_i1 = R::dnorm(A_draw(i,j),A_prior(i,j),tau0(i,j),false) * p_i;
          double u_i2 = R::dnorm(A_draw(i,j),A_prior(i,j),tau1(i,j),false) * (1-p_i);
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
      for(int pp=0; pp < (p+1); pp++){
        if(pp==0){
          if(exo){
            mat A_con = A_draw.rows(p*M, p*M+Mstar-1); 
            mat V_con = V_prior.rows(p*M, p*M+Mstar-1); 
            mat P_con = A_prior.rows(p*M, p*M+Mstar-1);
            int r = A_con.n_rows; int c = A_con.n_cols; int d = A_con.n_elem;
            
            lambda2_A(0,1) = sample_lambda2(V_con, A_tau(0,1), d_lambda, e_lambda, d, 1);
            sample_theta(V_con, A_con, P_con, lambda2_A(0,1), A_tau(0,1), r, c, false);
            V_prior.rows(p*M, p*M+Mstar-1) = V_con;
            
            if(sample_A){
              vec theta_vec = V_con.as_col(); vec lambda_vec = lambda2_A.submat(0,1,pp,1); double lambda_prod = prod(lambda_vec);
              sample_tau(A_tau(0,1), lambda_prod, theta_vec, A_tuning(0,1), A_accept(0,1), burnin, irep);
            }
          }
        }else{
          mat A_end = A_draw.rows((pp-1)*M, pp*M-1); 
          mat V_end = V_prior.rows((pp-1)*M, pp*M-1); 
          mat P_end = A_prior.rows((pp-1)*M, pp*M-1);
          int r_end = A_end.n_rows; int c_end = A_end.n_cols; int d_end = A_end.n_elem;
          
          double prodlambda = 1;
          if(pp>1){
            vec lambdavec  = lambda2_A.submat(1,0,pp-1,0);
            prodlambda = prod(lambdavec);
            }
          lambda2_A(pp,0) = sample_lambda2(V_end, A_tau(pp,0), d_lambda, e_lambda, d_end, prodlambda);
          sample_theta(V_end, A_end, P_end, lambda2_A(pp,0), A_tau(pp,0), r_end, c_end, false);
          V_prior.rows((pp-1)*M, pp*M-1) = V_end;
          
          if(sample_A){
            vec theta_vec_end = V_end.as_col();
            vec lambda_vec_end = lambda2_A.submat(1,0,pp,0); 
            double lambda_prod_end = prod(lambda_vec_end); 
            sample_tau(A_tau(pp,0), lambda_prod_end, theta_vec_end, A_tuning(pp,0), A_accept(pp,0), burnin, irep);
          }
          // weakly exogenous
          if(exo){
            mat A_exo = A_draw.rows(p*M+pp*Mstar, p*M+(pp+1)*Mstar-1); 
            mat V_exo = V_prior.rows(p*M+pp*Mstar, p*M+(pp+1)*Mstar-1); 
            mat P_exo = A_prior.rows(p*M+pp*Mstar, p*M+(pp+1)*Mstar-1);
            int r_exo = A_exo.n_rows; int c_exo = A_exo.n_cols; int d_exo = A_exo.n_elem;
            
            vec lambdavec = lambda2_A.submat(0,1,pp-1,1);
            double prodlambda = prod(lambdavec);
            lambda2_A(pp,1) = sample_lambda2(V_exo, A_tau(pp,1), d_lambda, e_lambda, d_exo, prodlambda);
            sample_theta(V_exo, A_exo, P_exo, lambda2_A(pp,1), A_tau(pp,1), r_exo, c_exo, false);
            V_prior.rows(p*M+pp*Mstar, p*M+(pp+1)*Mstar-1) = V_exo;
            
            if(sample_A){
              vec theta_vec_exo = V_exo.as_col();
              vec lambda_vec_exo = lambda2_A.submat(0,1,pp,1);
              double lambda_prod_exo = prod(lambda_vec_exo);
              sample_tau(A_tau(pp,1), lambda_prod_exo, theta_vec_exo, A_tuning(pp,1), A_accept(pp,1), burnin, irep);
            }
          }
        }
      }
      // coefficients H matrix
      int r = L_draw.n_rows; int c = L_draw.n_cols;
      lambda2_L = sample_lambda2(L_prior, L_tau, d_lambda, e_lambda, v, 1);
      sample_theta(L_prior, L_draw, l_prior, lambda2_L, L_tau, r, c, true);
      
      if(sample_A){
        vec theta_vec_l(v); int vv=0;
        for(int i=1; i < r; i++){
          for(int j=0; j < i; j++){
            theta_vec_l(vv) = L_prior(i,j); vv += 1;
          }
        }
        sample_tau(L_tau, lambda2_L, theta_vec_l, L_tuning, L_accept, burnin, irep);
      }
    }
    //-----------------------------------------------
    // Step 3: Sample covariances
    for(int mm=0; mm < M; mm++){
      vec data_sv = Em_str_draw.col(mm);
      vec cur_sv  = Sv_draw.unsafe_col(mm);  // changed to **unsafe**_col which reuses memory
      if(sv){
        const vec datastand = log(data_sv%data_sv + offset);
        double mu = Sv_para(0, mm);
        double phi = Sv_para(1, mm);
        double sigma = Sv_para(2, mm);
        double h0 = Sv_para(3, mm);
        stochvol::update_fast_sv(datastand, mu, phi, sigma, h0, cur_sv, rec, prior_spec, expert);
        Sv_para.col(mm) = arma::colvec({mu, phi, sigma});
        //Sv_draw.col(mm) = cur_sv;  // unsafe_col overwrites the original data without copying
      }else{
        sample_sig2(cur_sv, data_sv, a_1, b_1, T);
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
          //H_accmat(0,0) = H_accept;
          //accept_store.row((irep-burnin)/thin) = join_rows(A_accept,H_accmat);
          //H_tunmat(0,0) = H_tuning;
          //tuning_store.row((irep-burnin)/thin) = join_rows(A_tuning,H_tunmat);
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
