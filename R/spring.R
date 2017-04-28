##' Select Primary Relationships IN the General linear model.
##'
##' Adjust a multivariate regression with structuring, sparse penalty.
##'
##' @param x matrix of features. Do NOT include intercept.
##'
##' @param y matrix of responses.
##'
##' @param lambda1 sequence of decreasing \eqn{\ell_1}{l1}-penalty
##' levels. If \code{NULL} (the default), a vector is generated with
##' \code{nlambda1} entries, starting from a guessed level
##' \code{lambda1.max} where only the intercept is included, then
##' shrunken to \code{min.ratio*lambda1.max}.
##'
##' @param lambda2 real scalar; tunes the \eqn{\ell_2}{l2} structuring
##' penalty in the Elastic-net. Default is 0.05.
##'
##' @param struct matrix structuring the coefficients, possibly
##' sparsely encoded. If \code{NULL} (the default), the identity matrix is
##' used. See details below.
##'
##' @param intercept logical; indicates if a vector of intercepts
##' should be included in the model. Default is \code{TRUE}.
##'
##' @param normalize logical; indicates if predictor variables should
##' be normalized to have unit L2 norm before fitting.  Default is
##' \code{TRUE}.
##'
##' @param cov the matrix of variance-covariance between the
##' reponses y. If \code{NULL} (the default), will be inferred.
##'
##' @param nlambda1 integer that indicates the number of values to put
##' in the \code{lambda1} vector.  Ignored if \code{lambda1} is
##' provided.
##'
##' @param min.ratio minimal value of \eqn{\ell_1}{l1}-part of the
##' penalty that will be tried, as a fraction of the maximal
##' \code{lambda1} value. A too small value might lead to unstability
##' at the end of the solution path corresponding to small
##' \code{lambda1} combined with \eqn{\lambda_2=0}{lambda2=0}.  The
##' default value tries to avoid this, adapting to the
##' '\eqn{n<p}{n<p}' context. Ignored if \code{lambda1} is provided.
##'
##' @param verbose integer; activate verbose mode from '0' (nothing)
##' to '2' (detailed output).
##'
##' @param max.iter  integer; the maximal number of iteration
##' (i.e. number of alternated optimization between each parameter)
##' used to solve the problem for a given value of lambda1. Default is
##' 100
##'
##' @param threshold a threshold for convergence for each
##' \code{lambda1}. The algorithm stops the alternate scheme when the
##' norm of two successive values of the parameters the is below this
##' threshold. Default is \code{1e-3}.
##'
##' @param comp.df <...>
##'
##' @param mc.cores if \code{lambda2} is a vector, spring is run in
##' parallel: the default is to use as many core as there are entries
##' in \code{lambda2}, limited by the physical nulber of cores itself.
##'
##' @return an object with class \code{spring}, see the
##' documentation page \code{\linkS4class{spring}} for details.
##'
##' @seealso See also \code{\linkS4class{spring}} and \code{\link{cv.spring}}.
##' @name spring
##' @rdname spring
##' @keywords models, regression
##'
##' @export
spring <- function(x, y,
                   lambda1     = NULL , ## if NULL, an appropriate vector is generated
                   lambda2     = c(1,0.1,0.01,0.001),
                   struct      = sparseMatrix(i=1:p,j=1:p,x=rep(1,p)),
                   cov         = NULL , ## if NULL, will be inferred
                   intercept   = TRUE ,
                   normalize   = FALSE,
                   threshold   = 1e-3 ,
                   max.iter    = 100  ,
                   verbose     = ifelse(length(lambda2)>1,1,2) ,
                   min.ratio   = 1e-2 ,
                   nlambda1    = ifelse(is.null(lambda1),
                                 ifelse(is.null(cov),10,50),length(lambda1)),
                   comp.df     = TRUE,
                   mc.cores    = min(length(lambda2), detectCores())) {

  ## =============================================================
  ## INITIALIZATION & PARAMETERS RECOVERY
  if (Sys.info()[['sysname']] == "Windows") {
    warning("\nWindows does not support fork, enforcing mc.cores to '1'.")
    mc.cores <- 1
  }

  ## Initialization
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)

  ## not available when n < q
  stopifnot(n>q)

  ## Intercept treatment
  if (intercept) {
    xbar <- colMeans(x)
    ybar <- colMeans(y)
    x <- sweep(x, 2L, xbar, check.margin = FALSE)
    y <- sweep(y, 2L, ybar, check.margin = FALSE)
  }

  normy   <- sqrt(colSums(y^2))
  if (normalize) {
    normx <- sqrt(colSums(x^2)*n)
    x <- sweep(x, 2L, normx, "/", check.margin = FALSE)
#    y <- sweep(y, 2L, normy, "/", check.margin = FALSE)
  } else {
    normx <- rep(1,p)
#    normy <- rep(1,q)
  }

  ## Compute quantities that won't change along the path
  SXX    <- crossprod(x)/n
  SYY    <- crossprod(y)/n
  SYYm1  <- solve(SYY)
  SXY    <- crossprod(x,y)/n


  return(mclapply(lambda2, function(lbd2) {
    if (verbose>0) {cat("\n FITTING FOR LAMBDA2 =", lbd2, "\n")}
    out <- spring.learn(x=x, y=y, xbar=xbar,ybar,normx=normx,normy=normy,
                        SXX=SXX, SYY=SYY, SYYm1=SYYm1, SXY=SXY,
                        lambda1     = lambda1,
                        lambda2     = lbd2   ,
                        struct      = struct ,
                        cov         = cov    ,
                        intercept   = intercept,
                        normalize   = normalize,
                        threshold   = threshold,
                        max.iter    = max.iter ,
                        comp.df     = comp.df,
                        verbose     = verbose  ,
                        min.ratio   = min.ratio,
                        nlambda1    = nlambda1)
    return(out)
  }, mc.cores=mc.cores))

}

spring.learn <- function(x, y, xbar, ybar, normx, normy,
                         SXX       ,
                         SYY       ,
                         SYYm1     ,
                         SXY       ,
                         lambda1   ,
                         lambda2   ,
                         struct    ,
                         cov       ,
                         intercept ,
                         normalize ,
                         threshold ,
                         max.iter  ,
                         comp.df   ,
                         verbose   ,
                         monitoring,
                         min.ratio,
                         nlambda1) {


  ithreshold = 1e-3
  imax.iter  = 1000

  ## Compute quantities that won't change along the path
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)
  SXX.L2 <- SXX + lambda2 * struct

  ## Get an appropriate grid of l1-penalties, starting from the null model
  if (is.null(lambda1)) {
    lambda1 <- get.lambda1(SXY, nlambda1, min.ratio)
  } else {
    nlambda1 <- length(lambda1)
  }
  if (!inherits(struct, "Matrix"))
    struct <- as(struct, "dgCMatrix")

  ## check for oracle mode (i.e., when covariance of the response is known/provided)
  oracle <- ifelse (!(is.null(cov)), TRUE, FALSE)

  if (oracle) {
    ## ------------------------------------------------------------
    ## ORACLE MODE: FAST, WORKING WITH KNOW COVARIANCE
    ## Cost an Elastic-net/Lasso problem with (p x q) variables
    ##
    if (verbose >1)
      cat("\n\nCovariance of the output known: switching to 'Oracle' mode.")

    ## Call to learn.par function
    out <- learn.par.vec(SXX.L2, SXY, lambda1, cov, imax.iter, ithreshold, verbose)
    par.hat <- out$par.hat
    B.hat   <- out$B.hat
    cov.hat <- rep(list(cov),length(par.hat))
    inv.cov <- solve(cov)

    ## Post computation of the loglikelihood and the criterion optimized
    loglik <- loglikelihood(x,y,B.hat,inv.cov)
    pen.l1 <- sapply(par.hat, l1.norm)
    pen.l2 <- sapply(par.hat, l2.norm, struct, cov)
    crit   <- -1/n * loglik + lambda1 * pen.l1 + lambda2 * pen.l2 / 2
    monitor <- list(loglik = loglik, objective = crit)
  } else {
    ## ------------------------------------------------------------
    ## FULL INFERENCE: ALTERNATE OPTIMIZATION
    ##
    par.hat <- list()
    cov.hat <- list()
    B.hat   <- list()
    par.hat.c <- NULL
    monitor <- list()

    if (verbose >1) {
      cat("\n FIRST PASS: compute the initial covariance for the whole path")
    }
    cov.0 <- diag(diag(SYY))
    ## Compute the initial covariance estimator for the whole path of lambda1
    out <- learn.par.vec(SXX.L2, SXY, lambda1, diag(diag(SYY)), imax.iter, ithreshold, verbose)
    par.hat.0 <- out$par.hat
    cov.hat.0 <- lapply(par.hat.0, function(par) {
      return(learn.cov(SYY,SYYm1,SXX.L2,par))
    })

    if (verbose >1) {
      cat("\n")
      cat("\n Now run for something")
      cat("\n")
      cat("\n l1 penalty | objective")
    }
    for (l in 1:nlambda1) {

      lambda1.old <- ifelse(l>1,lambda1[l-1],lambda1[l])

      if (verbose >1) {cat("\n",format(lambda1[l],digits=3, width=10),"|")}

      ## Inference of the parameters
      cov.hat.c <- cov.hat.0[[l]]
      par.hat.c <- learn.par(SXX.L2, SXY, lambda1[l], lambda1.old, cov.hat.c$R, beta0=par.hat.0[[l]], maxit=imax.iter, thr=ithreshold)
      B.hat.c   <- -par.hat.c %*% cov.hat.c$R

      ## Evolution of the criterion
      loglik <- loglikelihood(x,y,B.hat.c,cov.hat.c$Oyy)
      pen.l1 <- l1.norm(par.hat.c)
      pen.l2 <- l2.norm(par.hat.c, struct, cov.hat.c$R)
      crit   <- -1/n * loglik + lambda1[l] * pen.l1 + lambda2 * pen.l2 / 2

      iter <- 0
      cond    <- FALSE
      if (verbose ){pb <- utils::txtProgressBar(char="+",style=1); }
      while (!cond) {
        iter <- iter + 1
        if (verbose >1) {
          utils::setTxtProgressBar(pb, value=iter/max.iter)
        }
        ## Inference of the parameters
        par.hat.o <- par.hat.c
        cov.hat.c <- learn.cov(SYY,SYYm1,SXX.L2,par.hat.c)
        par.hat.c <- learn.par(SXX.L2, SXY, lambda1[l], lambda1.old, cov.hat.c$R, beta0=par.hat.c, maxit=imax.iter, thr=ithreshold, B.hat=B.hat.c)
        B.hat.c   <- -par.hat.c %*% cov.hat.c$R

        ## Evolution of the criterion
        loglik.c <- loglikelihood(x,y,B.hat.c,cov.hat.c$Oyy)
        loglik   <- c(loglik,loglik.c)
        pen.l1   <- l1.norm(par.hat.c)
        pen.l2   <- l2.norm(par.hat.c, struct, cov.hat.c$R)
        crit     <- c(crit,-1/n * loglik.c + lambda1[l] * pen.l1 + lambda2 * pen.l2 / 2)

        delta.par <- sum((par.hat.c-par.hat.o)^2) / sum( (par.hat.c)^2)
        delta.par[is.nan(delta.par)] <- 0
        delta.crit <- (crit[iter+1] - crit[iter])
        if ((abs(delta.par) < threshold) || iter >= max.iter) {
          cond <- TRUE
        }
      }
      if (verbose >1) {
        close(pb)
        cat("            |",format(crit[iter+1],digits=3,width=9))
      }

      par.hat[[l]] <- par.hat.c
      cov.hat[[l]] <- symmpart(cov.hat.c$R)
      B.hat[[l]]   <- B.hat.c
      if (q > 1){
        colnames(par.hat[[l]]) <- colnames(B.hat[[l]]) <-
          colnames(cov.hat[[l]]) <- rownames(cov.hat[[l]]) <- colnames(y)
      }
      monitor[[l]] <- list(loglik = loglik, objective = crit)
    }
  }

  if (verbose >1) {
    cat("\n")
    cat("\n Post computation of some usual statistics")
  }
  ## Post computation of some usual statistics
  fitted    <- lapply(B.hat, function(beta) x %*% beta)
  residuals <- lapply(fitted , function(fit) fit-y)
  r.squared <- sapply(residuals, function(res) {
    return(1-colSums(res^2)/normy^2)
  })

  if (verbose >1) {
    cat("\n Degrees of freedom")
  }

  ## degrees of freedom (think about algebra to make this cheaper)
  if (comp.df) {
    df <- sapply(1:nlambda1, function(k) {
      A <- which(par.hat[[k]] != 0)
      if (length(A)>0) {
        ## From wikipedia on matrix products:
        ## A is mxn, B is pxq
        ## (A kron B )ij = A(1+floor((i-1/p)), 1+floor((j-1)/p)) * B(1+(i-1)mod p, 1 + (j-1) mod q)

        ## old way:
        ## M1 <- suppressMessages(kronecker(cov.hat[[k]],struct))[A,A]
        ## M2 <- suppressMessages(kronecker(cov.hat[[k]], SXX.L2))[A,A]
        ind.cov <- 1 + floor((A-1)/p)
        ind.pre <- 1 + (A-1) %% p
        M1 <- cov.hat[[k]][ind.cov,ind.cov] * struct[ind.pre,ind.pre]
        M2 <- cov.hat[[k]][ind.cov,ind.cov] * SXX.L2[ind.pre,ind.pre]
        return(length(A) - lambda2*trace(M1 %*% chol2inv(chol(symmpart(matrix(M2,length(A),length(A)))))))
      } else {
        return(0)
      }
    })
  } else {
    df <- as.numeric(rep(NA, nlambda1))
  }

  if (intercept) {
    mu <- lapply(B.hat, function(beta) ybar - crossprod(beta,xbar))
    df <- q + df
  } else {
    mu <- rep(list(Matrix(0,p,q)),nlambda1)
  }

  if (verbose >1) {
    cat("\n Normalize back to the original scale")
  }
  if (normalize) {
    B.hat <- lapply(B.hat, function(beta) {
      return(sweep(beta, 1L, normx, "/"))
    })
    par.hat <- lapply(par.hat, function(omega) {
      return(sweep(omega, 1L, normx, "/"))
    })
  }

  if (verbose >1) {
    cat("\n DONE.\n")
  }

  ## Then we are done
  return(new("spring",
             coef.direct = par.hat  ,
             coef.regres = B.hat    ,
             covariance  = cov.hat  ,
             intercept   = mu       ,
             lambda1     = lambda1  ,
             lambda2     = lambda2  ,
             df          = df       ,
             residuals   = residuals,
             fitted      = fitted   ,
             r.squared   = r.squared,
             normx       = normx    ,
             monitoring  = monitor  ))

}

learn.par.vec <- function(SXX.L2, SXY, lambda1.vec, cov, maxit=1000, thr=1e-4, verbose=verbose) {

  O.hat <- list()
  B.hat <- list()
  nlambda1 <- length(lambda1.vec)

  if (verbose >1) {
    cat("\n ")
    pb <- utils::txtProgressBar(min = 0, max = nlambda1, width=nlambda1, char=".")
    utils::setTxtProgressBar(pb, 1)
  }

  O.hat[[1]] <- learn.par(SXX.L2, SXY, lambda1.vec[1], lambda1.vec[1], cov, maxit=maxit, thr=thr)
  B.hat[[1]] <- -O.hat[[1]] %*% cov

  if (nlambda1 > 1) {
    for (k in 2:nlambda1) {
      if (verbose >1) {utils::setTxtProgressBar(pb, k)}
      O.hat[[k]] <- learn.par(SXX.L2, SXY, lambda1.vec[k], lambda1.vec[k-1], cov, beta0=O.hat[[k-1]], maxit=maxit, thr=thr, B.hat=B.hat[[k-1]])
      B.hat[[k]] <- -O.hat[[k]] %*% cov
    }
  }
  return(list(par.hat=O.hat,B.hat=B.hat))
}

learn.par <- function(xtx.L2, xty, lambda1, lambda1.old, cov, beta0=Matrix(0,p,q), relaxo=FALSE, maxit=1000, thr=1e-4, B.hat=NULL) {

  p <- ncol(xtx.L2)
  q <- ifelse(is.null(dim(cov)),1,ncol(cov))

  ## recursive strong rule to discard an entire column
  if (is.null(B.hat)) {
      B.hat <- - beta0 %*% cov
  }
  null <- apply( abs(xty - xtx.L2 %*% B.hat), 1, max) < (2 *lambda1-lambda1.old)
  ##    cat("\ndiscarded: ", sum(null)*q)
  ##    cat(" kept: ", sum(!null)*q)

  if (q == 1) {
    AtA <- as.matrix(cov * xtx.L2[!null, !null])
    Atb <- as.numeric(xty[!null])
  } else {
    AtA <- as.matrix(suppressMessages(kronecker(cov,xtx.L2[!null, !null])))
    Atb <- as.vector(xty[!null, ])
  }
  if (sum(!null) != 0) {
    beta.hat <- as.numeric(.Call("coordinate_lasso", as.vector(beta0[!null, ]), Atb, as.matrix(AtA), lambda1, thr, maxit, PACKAGE="spring")$xk)
    return(sparseMatrix(i = rep(which(!null),q), j = rep(1:q,each=sum(!null)), x = beta.hat, dims=c(p,q)))
  } else {
    return(beta0)
  }
  ## } else {
  ##   if (!is.null(beta0)) {
  ##     beta0 <- as.numeric(beta0)
  ##   } else {
  ##     beta0 <- rep(0,p*q)
  ##   }
  ##   if (is.null(dim(cov))) {
  ##     AtA <- as.matrix(cov * xtx.L2)
  ##     Atb <- as.numeric(xty)
  ##     q <- 1
  ##   } else {
  ##     AtA <- as.matrix(suppressMessages(kronecker(cov,xtx.L2)))
  ##     Atb <- as.vector(xty)
  ##     q <- ncol(cov)
  ##   }
  ##   beta.hat <- as.numeric(.Call("coordinate_lasso", beta0, Atb, AtA, lambda1, thr, maxit, PACKAGE="spring")$xk)
  ##   return(Matrix(beta.hat, p, q))
  ## }

  ## if (relaxo) {
  ##   A <- beta.hat != 0
  ##   if (sum(A)>0)
  ##     beta.hat[A] <- chol2inv(chol(AtA[A,A])) %*% Atb[A]
  ## }

}

learn.cov <- function(SYY, SYYm1, SXX, O.xy) {

  O.SXX.O <- crossprod(O.xy, SXX) %*% O.xy
  if (sum(abs(O.SXX.O)) < 1e-2) {
    R   <- SYY
    Oyy <- SYYm1
  }
  else if (ncol(SYY) == 1) {
    R <- as.numeric(0.5*(-1 + sqrt(1+4*SYY*O.SXX.O))/(O.SXX.O))
    Oyy <- 1/R
  }
  else {
    ## better with this
    eig <- eigen(nearPD(O.SXX.O, ensureSymmetry=TRUE)$mat %*% SYY)
    U   <- eig$vectors
    Um1 <- solve(U)
    xi  <- eig$values
    eta <- (1 + sqrt(1+4*xi))/2

    R   <- SYY %*% U %*% Diagonal(x=1/eta) %*% Um1
    Oyy <- U %*% Diagonal(x=eta) %*% Um1 %*% SYYm1
    ## R   <- SYY %*% U %*% Diagonal(x=1/eta) %*% Um1
    ## Oyy <- U %*% Diagonal(x=eta) %*% Um1 %*% SYYm1
  }

  return(list(R=R,Oyy=Oyy))
}

learn.default.args <- function(n,p,user) {

  if (is.null(user$lambda2)) {
    lambda2 <- c(1,0.1,0.01,0.001)
  } else {
    lambda2 <- user$lambda2
  }
  nlambda1 <- ifelse(is.null(user$lambda1),
                     ifelse(is.null(user$nlambda1),10,user$nlambda1)
                     ,length(user$lambda1))
  return(list(
    lambda1     = NULL ,
    lambda2     = lambda2,
    struct      = sparseMatrix(i=1:p,j=1:p,x=rep(1,p)),
    cov         = NULL ,
    intercept   = TRUE ,
    normalize   = FALSE,
    threshold   = 1e-3 ,
    max.iter    = 100  ,
    verbose     = 0    ,
    min.ratio   = 1e-2 ,
    nlambda1    = nlambda1,
    comp.df     = TRUE,
    mc.cores    = min(length(lambda2), detectCores())))
}

## learn.par.old <- function(x, y, lambda1, lambda2, cov, struct, beta0=NULL, warm.start=TRUE) {

##   n <- nrow(x) # number of input
##   p <- ncol(x) # number of input
##   q <- ncol(y) # number of output

##   ## Call to quadrupen with the appropriate tranformation of the data
##   if (q == 1) {
##     C <- sqrt(as.numeric(cov))
##     Cm1 <- 1/C
##     A <- x * C / sqrt(n)
##     b <- - y * Cm1  / sqrt(n)
##     struct <- struct * as.numeric(cov)
##   } else {
##     C   <- chol(cov)
##     Cm1 <- solve(C)
##     A <- kronecker(x, as(C,"matrix"))  / sqrt(n)
##     b <- - as.numeric(t( y %*% Cm1 ))  / sqrt(n)
##     if( inherits(struct, "sparseMatrix")) {
##       struct <- as(kronecker(struct,as(cov, "dgCMatrix")), "dgCMatrix")
##     } else {
##       struct <- kronecker(struct,cov)
##     }
##   }
##   if (!is.null(beta0) & warm.start & lambda2 > 1e-2) {
##     beta0 <- as.numeric(drop(beta0))
##   } else {
##     beta0 <- NULL
##   }

##   beta <- elastic.net(A, b,
##                       lambda1   = lambda1,
##                       lambda2   = lambda2,
##                       intercept = FALSE  , ## The following
##                       normalize = FALSE  , ## parameters
##                       naive     = TRUE   , ## should not be
##                       struct    = struct , ## change by the user
##                       max.feat  = ncol(A),
##                       beta0     = beta0  ,
##                       ## control = list(usechol  = warm.start,
##                       ## threshold = ifelse(warm.start,1e-3,1e-7)),
##                       checkargs = FALSE  )@coefficients

##   if (length(lambda1) == 1) {
##     return(Matrix(drop(beta), p, q, byrow=T))
##   } else {
##     return(apply(beta, 1, function(b) {
##       Matrix(drop(b), p, q, byrow=T)
##     }))
##   }

## }
