#define STRICT_R_HEADERS
#include <float.h>
#include <Rcpp.h>
#include <math.h>
#include <stdlib.h>
#include <GIGrvg.h>
using namespace Rcpp;

//' @name do_rgig
//' @noRd
//[[Rcpp::export]]
double do_rgig1(double lambda, double chi, double psi) {
  
  if (chi == 0){
    chi = DBL_MIN;
  }
  
  if ( !(R_FINITE(lambda) && R_FINITE(chi) && R_FINITE(psi)) ||
       (chi <  0. || psi < 0)      ||
       (chi == 0. && lambda <= 0.) ||
       (psi == 0. && lambda >= 0.) ) {
    throw std::bad_function_call();
  }
  
  double res;
  
  // circumvent GIGrvg in these cases
  if (chi < 11 * DBL_EPSILON) {
    /* special cases which are basically Gamma and Inverse Gamma distribution */
    if (lambda > 0.0) {
      res = R::rgamma(lambda, 2.0/psi);
    }
    else {
      res = 1.0/R::rgamma(-lambda, 2.0/chi); // fixed
    }
  }
  
  else if (psi < 11 * DBL_EPSILON) {
    /* special cases which are basically Gamma and Inverse Gamma distribution */
    if (lambda > 0.0) {
      res = R::rgamma(lambda, 2.0/psi);  // fixed
    }
    else {
      res = 1.0/R::rgamma(-lambda, 2.0/chi); // fixed
    }
    
  } else {
    SEXP (*fun)(int, double, double, double) = NULL;
    if (!fun) fun = (SEXP(*)(int, double, double, double)) R_GetCCallable("GIGrvg", "do_rgig");
    
    res = as<double>(fun(1, lambda, chi, psi));
  }
  
  return res;
}
