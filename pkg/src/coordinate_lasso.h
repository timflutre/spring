#ifndef _spring_COO_LASSO_H
#define _spring_COO_LASSO_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

// Include Armadillo / Rcpp / R to C/C++ basics
#include <string.h>
#include <sys/time.h>
#include <RcppArmadillo.h>

RcppExport SEXP coordinate_lasso(SEXP X0, SEXP XTY, SEXP XTX, SEXP PEN, SEXP THR, SEXP MAXIT) ;

#endif
