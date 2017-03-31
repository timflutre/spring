##' Cross-validation function for spring fitting methods.
##'
##' Function that computes K-fold cross-validated error of a
##' \code{spring} fit on a grid of \code{lambda1,lambda2}.
##'
##' @param x matrix of features, possibly sparsely encoded
##' (experimental). Do NOT include intercept.
##'
##' @param y response vector.
##'
##' @param K integer indicating the number of folds. Default is 10.
##'
##' @param folds list of \code{K} vectors that describes the folds to
##' use for the cross-validation. By default, the folds are randomly
##' sampled with the specified K. The same folds are used for each
##' values of \code{lambda2}.
##'
##' @param verbose logical; indicates if the progression should be
##' displayed. Default is \code{TRUE}.
##'
##' @param ... additional parameters to overwrite the defaults of the
##' spring fitting procedure. See the corresponding documentation
##' (\code{\link{spring}}.
##'
##' @return A list with two objects of class "criterion" and "model" for which 
##' \code{plot} methods are available.
##'
##' @seealso \code{\linkS4class{criterion}}, \code{\linkS4class{model}}, \code{\link{plot,criterion-method}}, \code{\link{plot,model-method}}.
##'
##'
##' @keywords models, regression
##' @name cv.spring
##' @aliases cv.spring
##' @rdname cv.spring
##'
##' @export
cv.spring <- function(x, y,
                      K        = 5,
                      folds    = split(sample(1:nrow(x)), rep(1:K, length=nrow(x))),
                      verbose  = TRUE,
                      ...) {

  n <- nrow(x)
  p <- ncol(x)
  user <- list(...)
  defs <- learn.default.args(n-max(sapply(folds,length)),p,user)
  args <- modifyList(defs, user)
  args$verbose <- FALSE
  args$comp.df <- FALSE ## save some time

  ## Compute a grid of lambda1 (the same for each fold)
  if (is.null(args$lambda1)) {
    args$lambda1 <- get.lambda1(crossprod(scale(x,args$intercept,args$normalize),scale(y,args$intercept,FALSE))/n, args$nlambda1,args$min.ratio)
  }

  ## CROSS-VALIDATION ON A GRID
  if (verbose){
    cat("\nCROSS-VALIDATION\n\n")
    cat(length(folds),"-fold CV on the lambda1 grid for each lambda2\n", sep="")
    cat("\tlambda2: ",args$lambda2)
    cat("\n\tomitting bloc")
  }

  nlambda2 <- length(args$lambda2)
  fold.err <- rep(list(matrix(NA, n, args$nlambda1)), nlambda2)

  ## CV work
  for (k in 1:K) {
    if (verbose){cat("", k)}
    omit <- folds[[k]]
    if(is.null(dim(y))) {
      y.train <- y[-omit]
      y.test <- y[omit]
    } else {
      y.train <- y[-omit, ]
      y.test <- y[omit, ]
    }

    ## Yep, this is a loop, but timing to save with vectorization is
    ## negligible compared to the time taken to fit the model
    fit <- do.call(spring, c(list(x=x[-omit, ], y=y.train), args))
    for(i in 1:nlambda2) {
      fold.err[[i]][omit,] <-
        sapply(predict(fit[[i]],  matrix(x[omit,], nrow=length(omit))), function(y.hat) {
          rowMeans((y.hat - y.test)^2, na.rm=TRUE)
        })
    }
  }

  mean <- lapply(fold.err, colMeans,  na.rm=TRUE)
  serr <- lapply(1:nlambda2, function(j) {
    sqrt((colSums(sweep(fold.err[[j]], 2L, mean[[j]], check.margin = FALSE)^2,na.rm=TRUE)/(n-1)) /K)
  })
  ## formatting cv.error for ggplot
  cv <- data.frame(value = unlist(mean), serr=unlist(serr),
                   lambda1 = rep(args$lambda1,nlambda2),
                   lambda2 = rep(args$lambda2, each=args$nlambda1))

  ## Recovering the best lambda1 and lambda2
  lambda1.cv <- sapply(1:nlambda2, function(j) {
    cv.min <- min(mean[[j]])
    lambda1.min   <- max(args$lambda1[mean[[j]] <= cv.min], na.rm=TRUE)
    lambda1.1se   <- max(args$lambda1[mean[[j]] <(mean[[j]]+serr[[j]]+1e-5)[match(lambda1.min,args$lambda1)]], na.rm=TRUE)
    return(c(cv.min, lambda1.min, lambda1.1se))
  })

  ind.lb2.min <- which.min(lambda1.cv[1,])
  lambda2     <- args$lambda2[ind.lb2.min]
  lambda1.min <- lambda1.cv[2, ind.lb2.min]
  lambda1.1se <- lambda1.cv[3, ind.lb2.min]

  ## Apply the fitting procedure with these best lambda2 parameter
  args$lambda2 <- lambda2
  if (verbose) {
    cat("\nFitting model with cross-validated parameters...\n\n")
  }
  best.fit <- do.call(spring, c(list(x=x,y=y),args))[[1]]

  ## Finally recover the CV choice (minimum and 1-se rule)
  ind.max <- length(best.fit@coef.direct)
  ind.min <- min(match(lambda1.min, args$lambda1),ind.max)
  ind.1se <- min(match(lambda1.1se, args$lambda1),ind.max)

  return(list(
      cv.min = list(
        criterion = new("criterion", name = "CV-error (min)", data = cv),
        model     = new("model",
                      lambda1      = lambda1.min,
                      lambda2      = lambda2,
                      coef.regres  = Matrix(best.fit@coef.regres[[ind.min]]),
                      coef.direct  = Matrix(best.fit@coef.direct[[ind.min]]),
                      covariance   = Matrix(best.fit@covariance [[ind.min]]),
                      intercept    = Matrix(best.fit@intercept  [[ind.min]]))),
      cv.1se = list(
        criterion = new("criterion", name = "CV-error (1se)", data = cv),
        model = new("model",
                      lambda1      = lambda1.1se,
                      lambda2      = lambda2    ,
                      coef.regres  = Matrix(best.fit@coef.regres[[ind.1se]]),
                      coef.direct  = Matrix(best.fit@coef.direct[[ind.1se]]),
                      covariance   = Matrix(best.fit@covariance [[ind.1se]]),
                      intercept    = Matrix(best.fit@intercept  [[ind.1se]])))))
}

