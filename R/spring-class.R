##' Class "spring"
##'
##' Class of object returned by the \code{spring} function.
##'
##' @slot coef.direct List of pxq Matrices (with class
##' \code{"dgCMatrix"}) of estimated direct effects in the model with
##' respect to the original input. The number of element in the list
##' corresponds the length of \code{lambda1}.
##'
##' @slot coef.regres List of pxq Matrices (with class
##' \code{"dgCMatrix"}) of estimated regression coefficients (both
##' direct and indirect effects) in the model with respect to the
##' original input. The number of element in the list corresponds the
##' length of \code{lambda1}.
##'
##' @slot covariance List of qxq matrices of estimated
##' covariance between the output. The number of element in the list
##' corresponds the length of \code{lambda1}.
##'
##' @slot intercept List of qx1 vectors of estimated intercept
##' terms. The number of element in the list corresponds the length of
##' \code{lambda1}.
##'
##' @slot lambda1 Vector (class \code{"numeric"}) of
##' \eqn{\ell_1}{l1} norm penalty levels for which the model has
##' eventually been fitted.
##'
##' @slot lambda2 Scalar for the level of (structuring)
##' \eqn{\ell_2}{l2} norm penalty applied.
##'
##' @slot df Estimated degree of freedoms for the successive
##' \code{lambda1}.
##'
##' @slot residuals List of Matrices of residuals. Each entry
##' of the list corresponds to a value of \code{lambda1}.}
##' \item{\code{r.squared}:}{Vector (class \code{"numeric"}) given the
##' coefficient of determination as a function of lambda1.
##'
##' @slot normx Vector (class \code{"numeric"}) containing the
##' square root of the sum of squares of each column of the design
##' matrix.
##'
##' @slot monitoring List (class \code{"list"}) which
##' contains evolution of the loglikelihood and of the objective
##' function for each lambda1 and each iteration of the algorithm.
##'
##' @slot fitted List of Matrices of fitted values. Each
##' entry of the list corresponds to a value of \code{lambda1}.
##'
##' @section Methods: This class comes with the usual
##' \code{predict(object, newx, ...)}, \code{fitted(object, ...)},
##' \code{residuals(object, ...)}, \code{print(object, ...)},
##' \code{show(object)} and \code{deviance(object, ...)} generic
##' (undocumented) methods.
##'
##' The two methods \code{loglik(object)} and \code{objective(object)}
##' are also defined to recover the log-likelihood and objective
##' functions along the fit. They can be applied on a spring object or
##' on a list of spring object. The corresponding output is an object
##' with class \code{criterion}. See the corrdesponding documentation
##' for defintiion and methods.
##'
##' The method \code{spr.model(object, lambda1, lambda2=NULL))}
##' applies on a spring object or a list of spring objects to extract
##' an single object with class \code{model}.  See the corresponding
##' documentation for definition and methods.
##'
##' @aliases fitted,spring-method predict,spring-method
##' deviance,spring-method print,spring-method
##' show,spring-method residuals,spring-method
##' loglik,spring-method objective,spring-method
##' objective,list-method loglik,list-method
##'
##' @docType class
##'
##' @keywords class
##'
##' @seealso See also \code{\link{spring}},
##' \code{\linkS4class{criterion}}, \code{\linkS4class{model}}.
##'
##' @name spring-class
##' @rdname spring-class
##'
##' @exportClass spring
##' @exportMethod fitted
##' @exportMethod residuals
##' @exportMethod predict
##' @exportMethod deviance
##' @exportMethod print
##' @exportMethod loglik
##' @exportMethod objective
##' @exportMethod show
##'
##' @importFrom stats fitted predict residuals deviance
setClassUnion("mu.cl", c("NULL","list"))
setClassUnion("matORnum", c("matrix", "numeric"))
setClass("spring",
  slots = c(
    coef.direct = "list"   ,
    coef.regres = "list"   ,
    covariance  = "list"   ,
    intercept   = "mu.cl"  ,
    lambda1     = "numeric",
    lambda2     = "numeric",
    df          = "numeric",
    residuals   = "list"   ,
    fitted      = "list"   ,
    r.squared   = "matORnum",
    normx       = "numeric",
    monitoring  = "list"   )
)

setMethod("print", "spring", definition =
   function(x, ...) {
     p <- nrow(x@coef.regres[[1]])
     q <- ncol(x@coef.regres[[1]])
     cat("Multivariate Linear regression with structured and sparse penalizer.\n")
     if (!is.null(x@intercept)) {
       cat("- number of predictors:", p,"+ intercept\n")
     } else {
       cat("- number of predictors:", p,"(no intercept)\n")
     }
       cat("- number of responses:", q,"\n")
     cat("- penalty parameter lambda1: ", length(x@lambda1), " points from ",
         format(max(x@lambda1), digits = 3)," to ",
         format(min(x@lambda1), digits = 3),"\n", sep="")
     cat("- penalty parameter lambda2: ", x@lambda2, "\n", sep="")

     invisible(x)
   }
)

setMethod("show", "spring", definition =
   function(object) {print(object)}
)

setMethod("fitted", "spring", definition =
   function(object, ...) {
     return(object@fitted)
   }
)

setMethod("predict", "spring", definition =
   function (object, newx=NULL, ...)  {
     if (is.null(newx)) {
       return(object@fitted)
     } else {
       if (is.null(object@intercept)) {
         return(lapply(object@coef.regres, function(coef) newx %*% coef))
       } else {
         return(
           lapply(1:length(object@coef.regres), function(k)
                  sweep(newx %*% object@coef.regres[[k]],2L,-object@intercept[[k]],check.margin=FALSE)))
       }

     }
   }
)

setMethod("residuals", "spring", definition =
   function(object, newx=NULL, newy=NULL, ...) {
     if (is.null(newx) | is.null(newy)) {
       return(object@residuals)
     } else {
       n <- length(object@lambda1)
       return(lapply(predict(object, newx), function(pred) newy - pred) )
     }
   }
)

setMethod("deviance", "spring", definition =
   function(object, newx=NULL, newy=NULL, ...) {
     dev <- sapply(residuals(object, newx, newy), function(res) {
       return(colSums(res^2))
     })
     return(dev)
   }
)

setGeneric("loglik", function(object) {standardGeneric("loglik")})
setMethod("loglik", "spring", definition =
   function(object) {
     if (is.null(object@monitoring$loglik)) { ## handle the oracle mode
       vloglik <- sapply(object@monitoring, function(x) x$loglik[length(x$loglik)])
     } else {
       vloglik <- object@monitoring$loglik
     }

     fraction <- sapply(object@coef.direct, function(coef) sum(abs(coef)))
     fraction <- fraction/max(fraction)

     return(new("criterion",
                name = "loglikelihood",
                data = data.frame(value = vloglik, lambda1 = object@lambda1, df=object@df,
                  fraction = fraction, lambda2 = rep(object@lambda2, length(object@lambda1)))))


   })

setMethod("loglik", "list", definition =
          function(object) {

            stopifnot(all(sapply(object, class)=="spring"))
            return(new("criterion", name = "loglikelihood",
                       data = do.call(rbind, lapply(object, function(obj) loglik(obj)@data))))

          })

setGeneric("objective", function(object) {standardGeneric("objective")})
setMethod("objective", "spring", definition =
   function(object) {
     if (is.null(object@monitoring$objective)) { ## handle the oracle mode
       vobjective <- sapply(object@monitoring, function(x) x$objective[length(x$objective)])
     } else {
       vobjective <- object@monitoring$objective
     }

     fraction <- sapply(object@coef.direct, function(coef) sum(abs(coef)))
     fraction <- fraction/max(fraction)

     return(new("criterion",
                name = "objective function",
                data = data.frame(value = vobjective, lambda1 = object@lambda1, df=object@df,
                  fraction = fraction, lambda2 = rep(object@lambda2, length(object@lambda1)))))

   })

setMethod("objective", "list", definition =
          function(object) {

            stopifnot(all(sapply(object, class)=="spring"))
            return(new("criterion", name = "objective function",
                       data = do.call(rbind, lapply(object, function(obj) objective(obj)@data))))

          })


##' Penalized criteria based on estimation of degrees of freedom
##'
##' Send back the values of some penalized criteria accompanied with
##' the vector(s) of parameters selected accordingly. The default
##' behavior recover AIC and BIC (with respective factor \eqn{2}{2}
##' and \eqn{\log(n)}{log(n)}) yet the user can specify any penalty.
##'
##' @usage pen.crit(object, penalty=setNames(c(2, log(N)), c("AIC","BIC")))
##'
##' @param object A list or a single object with class
##' \code{spring}. Typically, an output of the \code{spring} fitting
##' function.
##' @param penalty a vector with as many penalties a desired. The
##' default contains the penalty corresponding to the AIC and the BIC
##' (\eqn{2}{2} and \eqn{\log(n)}{log(n)}). Setting the "names"
##' attribute, as done in the default definition, leads to outputs
##' which are easier to read.
##' @return A list of object with class \code{criterion} for which a
##' plot is available.
##' @seealso \code{\linkS4class{spring}}.
##'
##' @references Ryan Tibshirani and Jonathan Taylor. Degrees of
##' freedom in lasso problems, Annals of Statistics, 40(2) 2012.
##'
##' @name pen.crit
##' @aliases pen.crit pencrit,spring-method
##' @docType methods
##' @rdname pen.crit
##'
##' @import ggplot2 reshape2 scales grid methods
##' @exportMethod pen.crit
setGeneric("pen.crit", function(object, penalty=stats::setNames(c(2, log(N)), c("AIC","BIC"))) {
             N <- unique(sapply(object, function(obj) nrow(residuals(obj)[[1]])))
             standardGeneric("pen.crit")
           })

setMethod("pen.crit", "spring", definition =
          function(object, penalty=stats::setNames(c(2, log(N)), c("AIC","BIC"))) {
            N <- unique(sapply(object, function(obj) nrow(residuals(obj)[[1]])))
            return(pen.crit(list(object)))
          })

setMethod("pen.crit", "list", definition =
          function(object, penalty=stats::setNames(c(2, log(N)), c("AIC","BIC"))) {

            ## Check that all objects passed to the function are spring objects
            stopifnot(all(sapply(object, class)=="spring"))
            N <- unique(sapply(object, function(obj) nrow(residuals(obj)[[1]])))
            if (length(N) != 1) {stop("\nThe same data sould be use for every single spring fits.")}

            return(stats::setNames(lapply(1:length(penalty),
                          function(ipen) {
                            get.criteria(penalty[ipen], object, names(penalty)[ipen])
                          }
                          ), names(penalty)))
          })

get.one.criterion <- function(object, penalty) {

  lambda1 <- object@lambda1
  ## compute the penalized criteria

  value <- -2*loglik(object)@data$value + penalty * object@df

  fraction <- sapply(object@coef.direct, function(coef) sum(abs(coef)))
    fraction <- fraction/max(fraction)

  ## put together all relevant information about those criteria
  return(data.frame(value, df=object@df, lambda1=lambda1, fraction = fraction,
                    lambda2=rep(object@lambda2,length(lambda1))))
}

get.criteria <- function(penalty, objects, name) {

  crit <- lapply(objects, get.one.criterion, penalty)
  best <- which.min(lapply(crit, function(l) min(l$value) ))
  imin <- which.min(crit[[best]]$value)

  return(list(
      criterion = new("criterion",
          name        = name,
          data        = as.data.frame(do.call(rbind, crit))),
      model = model$new(
      lambda1     = objects[[best]]@lambda1[imin],
      lambda2     = objects[[best]]@lambda2      ,
      coef.regres = Matrix(objects[[best]]@coef.regres[[imin]]),
      coef.direct = Matrix(objects[[best]]@coef.direct[[imin]]),
      covariance  = Matrix(objects[[best]]@covariance [[imin]]),
      intercept   = Matrix(objects[[best]]@intercept  [[imin]]))))
}

##' Picking a model among a spring fit
##'
##' Send a model corresponding to a given couple of \code{lambda1},
##' \code{lambda2}.  If applied on a spring object, \code{lambda2} is
##' ignored. if applied on a list of spring object, sleection is made
##' according both \code{lambda1} and \code{lambda2}.
##'
##' @usage model(object, lambda1, lambda2=NULL)
##'
##' @param object A list or a single object with class
##' \code{spring}. Typically, an output of the \code{spring} fitting
##' function.
##' @param lambda1 a scalar. The model with the closest \code{lambda1} is picked out from the available fit in \code{object}.
##' @param lambda2 a scalar. The model with the closest \code{lambda2}
##' is picked out from the available fit in \code{object}. Useless
##' (and ignored) if a single \code{spring} object is provided.
##'
##' @return A object with class \code{model}. See the corresponding
##' documentation and method.
##' @seealso \code{\linkS4class{model}}.
##'
##' @name model
##' @aliases model model.spring model,spring-method model,list-method
##' @docType methods
##' @rdname model
##'
##' @exportMethod getModel
setGeneric("getModel", function(object, lambda1, lambda2=NULL)
           {standardGeneric("getModel")})

setMethod("getModel", "spring", definition =
          function(object, lambda1, lambda2=NULL) {

            ilambda1 <- which.min(abs(lambda1 - object@lambda1))

            return(model$new(
                       lambda1      = object@lambda1[[ilambda1]],
                       lambda2      = object@lambda2,
                       coef.regres  = Matrix(object@coef.regres[[ilambda1]]),
                       coef.direct  = Matrix(object@coef.direct[[ilambda1]]),
                       covariance   = Matrix(object@covariance [[ilambda1]]),
                       intercept    = Matrix(object@intercept  [[ilambda1]]))
                   )

          })

setMethod("getModel", "list", definition =
          function(object, lambda1, lambda2=NULL) {

            ## Check that all objects passed to the function are spring objects
            stopifnot(all(sapply(object, class)=="spring"))

            if(is.null(lambda2)){stop("Lambda2 must be specified when a list of spring object is given")}

            ilambda2 <- which.min(abs(lambda2 - sapply(object, function(obj) obj@lambda2)))
            print(ilambda2)
            return(getModel(object[[ilambda2]], lambda1))

          })
