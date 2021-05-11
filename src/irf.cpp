// [[Rcpp::depends(RcppParallel,RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "helper.h"
using namespace Rcpp;
using namespace arma;

struct IrfParallel : public RcppParallel::Worker
{
  // inputs
  const cube A_large, S_large, Ginv_large;
  const List shocklist;
  const int type, nhor;
  Function Rchol;
  
  // destination matrix
  RcppParallel::RMatrix<double> irf_output;
  RcppParallel::RMatrix<double> rot_output;
  RcppParallel::RMatrix<int> counter_output;
  
  // initialize with source and destination
  IrfParallel(cube A_large, cube S_large, cube Ginv_large, List shocklist, int type, int nhor, Function Rchol, 
              NumericMatrix irf_output, NumericMatrix rot_output, IntegerMatrix counter_output) 
    : A_large(A_large), S_large(S_large), Ginv_large(Ginv_large), shocklist(shocklist), type(type), nhor(nhor), 
      Rchol(Rchol), irf_output(irf_output), rot_output(rot_output), counter_output(counter_output) {}
  
  // take the square root of the range of elements requested
  void operator()(std::size_t begin, std::size_t end) {
    for(int irep = begin; irep < end; irep++){
      // get stuff
      List shock_idx = shocklist["shock.idx"];
      LogicalVector shock_cidx = shocklist["shock.cidx"];
      const int plag = shocklist["plag"];
      // get parameter
      const int bigK = A_large.n_cols;
      const int N = shock_idx.size();
      // current draw
      mat Amat = A_large.row(irep);
      mat Smat = S_large.row(irep);
      mat Ginv = Ginv_large.row(irep);
      // construct Fmat
      cube Fmat(bigK, bigK, plag, fill::zeros);
      for(uword pp = 0; pp < plag; pp++){
        Fmat.slice(pp) = Amat.cols(pp*bigK,(pp+1)*bigK-1);
      }
      // create P0G
      mat P0G(bigK,bigK, fill::zeros);
      bool chol_success = TRUE;
      for(uword cc = 0; cc < N; cc++){
        uvec idx = shock_idx[cc];
        if(shock_cidx[cc]){
          // fall back on pivoting if Cholesky do not work
          mat Sig_chol; mat sigma = Smat.submat(idx,idx);
          chol_success = chol(Sig_chol, sigma, "lower");
          if(!chol_success){
            break;
          }
          /*
           // Fall back on Rs chol if armadillo fails (it suppports pivoting)
           if(chol_success == false){
           NumericMatrix tmp = Rchol(sigma, true, false, -1);
           int d = sigma.n_cols;
           mat cholV_tmp = mat(tmp.begin(), d, d, false);
           uvec piv = sort_index(as<vec>(tmp.attr("pivot")));
           Sig_chol = cholV_tmp.cols(piv);
           Sig_chol = Sig_chol.t();
           }
           */
          P0G.submat(idx,idx) = Sig_chol;
        }else{
          P0G.submat(idx,idx) = Smat.submat(idx,idx);
        }
      }
      if(chol_success){
        // create PHI
        cube PHI(bigK, bigK, nhor); get_PHI(PHI, Fmat, nhor);
        // find rotation matrix
        mat Q_bar = eye(bigK,bigK);
        if(type==3){
          int MaxTries = shocklist["MaxTries"];
          mat S_cube = shocklist["S.cube"];
          mat P_cube = shocklist["P.cube"];
          cube Z_cube = shocklist["Z.cube"];
          uvec shock_order = shocklist["shock.order"];
          uvec shock_horz = shocklist["shock.horz"];
          LogicalVector nozero = shocklist["no.zero.restr"];
          int H_restr = shock_horz.n_elem;
          int N_restr = bigK * H_restr;
          // build irf_restr
          mat irf_restr(0, bigK);
          for(uword hh = 0; hh < H_restr; hh++){
            uword horz = shock_horz(hh);
            mat irf_hh = PHI.slice(horz) * Ginv * P0G;
            irf_restr = join_cols(irf_restr, irf_hh);
          }
          // sort Zcube
          cube Z_cube_sorted(size(Z_cube), fill::zeros);
          for(uword i = 0; i < bigK; i++){
            uword idx = shock_order(i);
            Z_cube_sorted.slice(i) = Z_cube.slice(idx);
          }
          // set up while loop
          int icounter = 0;
          double condall = 0.0;
          while( condall == 0.0 && icounter < MaxTries){
            mat Q(bigK, bigK, fill::eye);
            for(uword cc = 0; cc < N; cc++){
              uvec idx = shock_idx[cc];
              int Kidx = idx.size();
              if(shock_cidx[cc]){
                mat randMat(Kidx, Kidx, fill::randn), Qc(Kidx, Kidx, fill::zeros);
                if(nozero[cc]){
                  mat R; qr(Qc, R, randMat);
                }else{
                  for(uword i = 0; i < Kidx; i++){
                    mat Ztemp = Z_cube_sorted.slice(idx[i]);
                    colvec Zsum = sum(abs(Ztemp),1); uvec zidx = find(Zsum);
                    Ztemp = Ztemp.rows(zidx);
                    mat R(0,Kidx);
                    if(i == 0){
                      mat R1 = Ztemp * irf_restr.cols(idx);
                      R = join_cols(R,R1);
                    }else{
                      mat R2 = Qc.head_cols(i).t();
                      if(Ztemp.n_elem != 0){
                        mat R1 = Ztemp * irf_restr.cols(idx);
                        R = join_cols(R, R1);
                      }
                      R = join_cols(R, R2);
                    }
                    mat Rt = R.t();
                    mat NU; get_nullspace(NU, Rt); 
                    vec x_j = randMat.col(i);
                    double div = arma::as_scalar((NU.t() * x_j).t() * (NU.t() * x_j));
                    vec q_j = NU * ( NU.t() * x_j / sqrt(div));
                    Qc.col(i) = q_j;
                  }
                }
                Q.submat(idx,idx) = Qc;
              }
            }
            // reorder again
            for(uword i = 0; i < bigK; i++){
              uword idx = shock_order(i);
              Q_bar.col(idx) = Q.col(i);
            }
            // create irfcheck
            mat irf_check = irf_restr * Q_bar;
            vec signCheck(bigK,1);
            for(uword kk = 0; kk < bigK; kk++){
              vec STemp = S_cube.col(kk);
              vec PTemp = P_cube.col(kk);
              if(sum(abs(STemp))>0){
                vec prob(N_restr); 
                for(uword nn = 0; nn < N_restr; nn++){
                  if(PTemp(nn) > R::runif(0,1)) prob(nn) = 1; else prob(nn) = 0;
                }
                mat PDiag = diagmat(prob);
                vec IrfTemp = sign(irf_check.col(kk));
                double getsum = arma::as_scalar(IrfTemp.t() * PDiag * STemp);
                double chksum = arma::as_scalar(sum(abs(PDiag * STemp)));
                if(getsum == chksum) signCheck(kk) = 1; else signCheck(kk) = 0;
              }else{
                signCheck(kk) = 1;
              }
            }
            condall = prod(signCheck);
            icounter += 1;
          }
          // save counter
          counter_output(irep,0) = icounter;
        }else{
          counter_output(irep,0) = 1;
        }
        // compute shock
        mat invGSigma_u = Ginv * P0G * Q_bar;
        // compute impulse responses
        cube irfa(bigK, bigK, nhor);
        for(uword ihor = 0; ihor < nhor; ihor++){
          irfa.slice(ihor) = PHI.slice(ihor) * invGSigma_u;
        }
        // transform to matrix
        const size_t chunksize = bigK*nhor;
        for(size_t ihor = 0; ihor < nhor; ihor++){
          const size_t ihor_begin = irep*chunksize + bigK*ihor;
          const size_t rot_begin = irep*bigK;
          
          for(size_t i = 0; i < bigK; i++){
            for(size_t j = 0; j < bigK; j++){
              irf_output(ihor_begin + i,j) = irfa(i,j,ihor);
              rot_output(rot_begin + i, j) = Q_bar(i,j);
            }
          }
        }
        // end of impulse response computation
      } // end choll_success
      // check user interruption
      Rcpp::checkUserInterrupt();
    }
  }
};

//' @name compute_irf_parallel
//' @noRd
//' @export
// [[Rcpp::export]]
List compute_irf_parallel(arma::cube A_large, arma::cube S_large, arma::cube Ginv_large, const int type, const int nhor, const int thindraws, const SEXP shocklist_in) {
  
  // input stuff
  List shocklist(clone(shocklist_in));
  const int bigK = A_large.n_cols;
  
  // Import Rs chol function
  Environment base = Environment("package:base");
  Function Rchol = base["chol"];
  
  // allocate the output matrix
  NumericMatrix irf_output(bigK*thindraws*nhor, bigK);
  std::fill(irf_output.begin(), irf_output.end(), NumericVector::get_na());
  NumericMatrix rot_output(bigK*thindraws, bigK);
  IntegerMatrix counter_output(thindraws,1);
  // set up worker
  IrfParallel irf(A_large, S_large, Ginv_large, shocklist, type, nhor, Rchol, irf_output, rot_output, counter_output);
  // call parallelFor to do the work
  RcppParallel::parallelFor(0, thindraws, irf, 100);
  
  // transform to cube for output
  cube irf_out(bigK*nhor, bigK, thindraws, fill::zeros);
  cube rot_out(bigK, bigK, thindraws, fill::zeros);
  const size_t chunksize = nhor*bigK;
  for(uword irep = 0; irep < thindraws; irep++){
    const size_t irf_begin = irep*chunksize;
    const size_t rot_begin = irep*bigK;
    
    for(size_t j = 0; j < bigK; j++){
      for(size_t i = 0; i < chunksize; i++){
        irf_out(i,j,irep) = irf_output(irf_begin + i, j);
      }
      for(size_t i = 0; i < bigK; i++){
        rot_out(i,j,irep) = rot_output(rot_begin + j, i);
      }
    }
  }
  // return the output matrix
  return List::create(Named("irf") = irf_out,
                      Named("rot") = rot_out,
                      Named("counter") = counter_output);
}

//' @name compute_irf
//' @noRd
// [[Rcpp::export]]
List compute_irf(arma::cube A_large, arma::cube S_large, arma::cube Ginv_large, const int type, const int nhor, const int thindraws, const SEXP shocklist_in) {
  
  // input stuff
  List shocklist(clone(shocklist_in));
  //----------------------------------------------------------------------------
  // miscellaneous
  List shock_idx = shocklist["shock.idx"];
  LogicalVector shock_cidx = shocklist["shock.cidx"];
  // parameter
  const int plag = shocklist["plag"];
  const int bigK = A_large.n_cols;
  const int N = shock_idx.size();
  const size_t chunksize = nhor*bigK;
  
  // Import Rs chol function
  Environment base = Environment("package:base");
  Function Rchol = base["chol"];
  
  // allocate the output matrix
  NumericMatrix irf_output(bigK*thindraws*nhor, bigK); 
  std::fill(irf_output.begin(), irf_output.end(), NumericVector::get_na());
  NumericMatrix rot_output(bigK*thindraws, bigK);
  IntegerMatrix counter_output(thindraws,1);
  //----------------------------------------------------------------------------
  for(int irep = 0; irep < thindraws; irep++){
    // current draw
    mat Amat = A_large.row(irep);
    mat Smat = S_large.row(irep);
    mat Ginv = Ginv_large.row(irep);
    // construct Fmat
    cube Fmat(bigK, bigK, plag, fill::zeros);
    for(uword pp = 0; pp < plag; pp++){
      Fmat.slice(pp) = Amat.cols(pp*bigK,(pp+1)*bigK-1);
    }
    // create P0G
    bool chol_success = TRUE;
    mat P0G(bigK,bigK, fill::zeros);
    for(uword cc = 0; cc < N; cc++){
      uvec idx = shock_idx[cc];
      if(shock_cidx[cc]){
        // fall back on pivoting if Cholesky do not work
        mat Sig_chol; mat sigma = Smat.submat(idx,idx);
        chol_success = chol(Sig_chol, sigma, "lower");
        if(!chol_success){
          break;
        }
        /*
        // Fall back on Rs chol if armadillo fails (it suppports pivoting)
        if(chol_success == false){
          NumericMatrix tmp = Rchol(sigma, true, false, -1);
          int d = sigma.n_cols;
          mat cholV_tmp = mat(tmp.begin(), d, d, false);
          uvec piv = sort_index(as<vec>(tmp.attr("pivot")));
          Sig_chol = cholV_tmp.cols(piv);
          Sig_chol = Sig_chol.t();
        }
         */
        P0G.submat(idx,idx) = Sig_chol;
      }else{
        P0G.submat(idx,idx) = Smat.submat(idx,idx);
      }
    }
    // iterate through if cholesky fails
    if(!chol_success){
      continue;
    }
    
    // create PHI
    cube PHI(bigK, bigK, nhor); get_PHI(PHI, Fmat, nhor);
    // find rotation matrix
    mat Q_bar = eye(bigK,bigK);
    if(type==3){
      int MaxTries = shocklist["MaxTries"];
      mat S_cube = shocklist["S.cube"];
      mat P_cube = shocklist["P.cube"];
      cube Z_cube = shocklist["Z.cube"];
      uvec shock_order = shocklist["shock.order"];
      uvec shock_horz = shocklist["shock.horz"];
      LogicalVector nozero = shocklist["no.zero.restr"];
      int H_restr = shock_horz.n_elem;
      int N_restr = bigK * H_restr;
      // build irf_restr
      mat irf_restr(0, bigK);
      for(uword hh = 0; hh < H_restr; hh++){
        uword horz = shock_horz(hh);
        mat irf_hh = PHI.slice(horz) * Ginv * P0G;
        irf_restr = join_cols(irf_restr, irf_hh);
      }
      // sort Zcube
      cube Z_cube_sorted(size(Z_cube), fill::zeros);
      for(uword i = 0; i < bigK; i++){
        uword idx = shock_order(i);
        Z_cube_sorted.slice(i) = Z_cube.slice(idx);
      }
      // set up while loop
      int icounter = 0;
      double condall = 0.0;
      while( condall == 0.0 && icounter < MaxTries){
        mat Q(bigK, bigK, fill::eye);
        for(uword cc = 0; cc < N; cc++){
          uvec idx = shock_idx[cc];
          int Kidx = idx.size();
          if(shock_cidx[cc]){
            mat randMat(Kidx, Kidx, fill::randn), Qc(Kidx, Kidx, fill::zeros);
            if(nozero[cc]){
              mat R; qr(Qc, R, randMat);
            }else{
              for(uword i = 0; i < Kidx; i++){
                mat Ztemp = Z_cube_sorted.slice(idx[i]);
                colvec Zsum = sum(abs(Ztemp),1); uvec zidx = find(Zsum);
                Ztemp = Ztemp.rows(zidx);
                mat R(0,Kidx);
                if(i == 0){
                  mat R1 = Ztemp * irf_restr.cols(idx);
                  R = join_cols(R,R1);
                }else{
                  mat R2 = Qc.head_cols(i).t();
                  if(Ztemp.n_elem != 0){
                    mat R1 = Ztemp * irf_restr.cols(idx);
                    R = join_cols(R, R1);
                  }
                  R = join_cols(R, R2);
                }
                mat Rt = R.t();
                mat NU; get_nullspace(NU, Rt); 
                vec x_j = randMat.col(i);
                double div = arma::as_scalar((NU.t() * x_j).t() * (NU.t() * x_j));
                vec q_j = NU * ( NU.t() * x_j / sqrt(div));
                Qc.col(i) = q_j;
              } // end inner for-loop
            } // end inner if-cond
            Q.submat(idx,idx) = Qc;
          } // end outer if-cond
        } // end outer for-loop
        // reorder again
        for(uword i = 0; i < bigK; i++){
          uword idx = shock_order(i);
          Q_bar.col(idx) = Q.col(i);
        }
        // create irfcheck
        mat irf_check = irf_restr * Q_bar;
        vec signCheck(bigK,1);
        for(uword kk = 0; kk < bigK; kk++){
          vec STemp = S_cube.col(kk);
          vec PTemp = P_cube.col(kk);
          if(sum(abs(STemp))>0){
            vec prob(N_restr);
            for(uword nn = 0; nn < N_restr; nn++){
              if(PTemp(nn) > R::runif(0,1)) prob(nn) = 1; else prob(nn) = 0;
            }
            mat PDiag = diagmat(prob);
            vec IrfTemp = sign(irf_check.col(kk));
            double getsum = arma::as_scalar(IrfTemp.t() * PDiag * STemp);
            double chksum = arma::as_scalar(sum(abs(PDiag * STemp)));
            if(getsum == chksum) signCheck(kk) = 1; else signCheck(kk) = 0;
          }else{
            signCheck(kk) = 1;
          }
        }
        condall = prod(signCheck);
        icounter += 1;
      }
      // save counter
      counter_output(irep,0) = icounter;
    }else{
      counter_output(irep,0) = 1;
    } // end if-cond type==3
    
    // compute shock
    mat invGSigma_u = Ginv * P0G * Q_bar;
    
    // compute impulse responses
    cube irfa(bigK, bigK, nhor);
    for(uword ihor = 0; ihor < nhor; ihor++){
      irfa.slice(ihor) = PHI.slice(ihor) * invGSigma_u;
    }
    // transform to matrix
    for(size_t ihor = 0; ihor < nhor; ihor++){
      const size_t ihor_begin = irep*chunksize + bigK*ihor;
      const size_t rot_begin = irep*bigK;
      
      for(size_t i = 0; i < bigK; i++){
        for(size_t j = 0; j < bigK; j++){
          irf_output(ihor_begin + i,j) = irfa(i,j,ihor);
          rot_output(rot_begin + i, j) = Q_bar(i,j);
        }
      }
    }
    // check user interruption
    if(irep % 200 == 0)
      Rcpp::checkUserInterrupt();
  } // end for-loop of impulse response computation
  //----------------------------------------------------------------------------
  // transform to cube for output
  cube irf_out(bigK*nhor, bigK, thindraws, fill::zeros);
  cube rot_out(bigK, bigK, thindraws, fill::zeros);
  for(size_t irep = 0; irep < thindraws; irep++){
    size_t irf_begin = irep*chunksize;
    size_t rot_begin = irep*bigK;
    
    for(size_t j = 0; j < bigK; j++){
      for(size_t i = 0; i < chunksize; i++){
        irf_out(i,j,irep) = irf_output(irf_begin + i, j);
      }
      for(size_t i = 0; i < bigK; i++){
        rot_out(i,j,irep) = rot_output(rot_begin + j, i);
      }
    }
  }
  // return the output matrix
  return List::create(Named("irf") = irf_out,
                      Named("rot") = rot_out,
                      Named("counter") = counter_output);
}
