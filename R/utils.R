## ' Generate a matrix for the discrete derivative operator
##
##' Taken from genlasso package. Usefull to construct the structring
##' matrix
##
##' @param p an integer for the number of features
##'
##' @param K an integer for the order (default is 1) for the order of
##' the derivative operator
##'
##' @return an sparse matrix with class \code{dgCMatrix}.
##'
##' @name diffOp
##' @rdname diffOp
##'
##' @export
diffOp <- function (p, K=1) {
  D0  <- bandSparse(p, p, 0:1, diagonals = list(rep(-1, p), rep(1, p - 1)))
  if (K <= 1) {
    D <- D0
  } else {
    D <- D0
    for (k in 1:(K-1)) {
      D <- D %*% D0
    }
  }
  return(D[1:(p-K), ])
}

get.lambda1 <- function(xty,nlambda1,min.ratio) {
  lmax <- max(abs(xty))
  return(10^seq(log10(lmax), log10(min.ratio*lmax), len=nlambda1))
}

trace <- function(A) {
  return(sum(diag(A)))
}

l1.norm <- function(x) {
  return(sum(abs(x)))
}

l2.norm <- function(x,L,R) {
  ## trace(crossprod(x,L) %*% tcrossprod(x,R))
  return(crossprod( as.vector(crossprod(x,L)),as.vector(tcrossprod(x,R)) ))
}

loglikelihood <- function(X,Y,B,Oyy) {

  n <- nrow(X)
  ldO <- log(det(Oyy))
  loglik <- function(Beta) {
    term <- as.matrix(Y-X %*% Beta)
    return(n/2 * ldO - 1/2 * trace(Oyy %*% crossprod(term) ))
  }
  ## n/2 * ldO - 1/2 trace(Oyy %*% (SYY - 2 SYX %*% Beta + t(Beta) %*% SXX %*% Beta)  )
  if (is.list(B)) {
    return(sapply(B, loglik))
  } else {
    return(loglik(B))
  }
}


