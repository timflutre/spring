#include "coordinate_lasso.h"

using namespace Rcpp;
using namespace arma;

SEXP coordinate_lasso(SEXP X0, SEXP XTY, SEXP XTX, SEXP PEN, SEXP THR, SEXP MAXIT) {
  
  vec x0   = as<vec>(X0)  ;
  mat xtx  = as<mat>(XTX) ;
  vec xty  = as<vec>(XTY) ;
  double pen = as<double>(PEN);
  double thr = as<double>(THR) ;
  vec diag = xtx.diag()   ;
  int max_iter = as<uword>(MAXIT) ; // max. number of iteration

  vec xtxw = xtx * x0;
  
  vec xk = x0          ; // output vector
  int j, i = 0         ; // current iterate
  double delta = 2*thr ; // change in beta
  double u, d          ; // temporary scalar
  
  while ((delta > thr) && (i < max_iter)) {
    delta = 0;
    for (j=0; j<x0.n_elem; j++) {
      // Soft thresholding operator
      u = x0(j) * diag(j) - xty(j) - xtxw(j);
      xk(j)  = fmax(1-pen/fabs(u),0) * u/diag(j) ;
      d = xk(j)-x0(j);
      delta += pow(d,2);
      xtxw  += d*xtx.col(j) ;
    }

    // preparing next iterate
    delta = sqrt(delta);
    x0 = xk;
    i++;
    
    R_CheckUserInterrupt();
  }

  return List::create(Named("xk")   = xk,
		      Named("xtxw") = xtxw);
  
}
