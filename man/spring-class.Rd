\docType{class}
\name{spring-class}
\alias{deviance,spring-method}
\alias{fitted,spring-method}
\alias{loglik,list-method}
\alias{loglik,spring-method}
\alias{objective,list-method}
\alias{objective,spring-method}
\alias{predict,spring-method}
\alias{print,spring-method}
\alias{residuals,spring-method}
\alias{show,spring-method}
\alias{spring-class}
\title{Class "spring"}
\description{
Class of object returned by the \code{spring} function.
}
\section{Slots}{

\describe{
\item{\code{coef.direct}}{List of pxq Matrices (with class
\code{"dgCMatrix"}) of estimated direct effects in the model with
respect to the original input. The number of element in the list
corresponds the length of \code{lambda1}.}

\item{\code{coef.regres}}{List of pxq Matrices (with class
\code{"dgCMatrix"}) of estimated regression coefficients (both
direct and indirect effects) in the model with respect to the
original input. The number of element in the list corresponds the
length of \code{lambda1}.}

\item{\code{covariance}}{List of qxq matrices of estimated
covariance between the output. The number of element in the list
corresponds the length of \code{lambda1}.}

\item{\code{intercept}}{List of qx1 vectors of estimated intercept
terms. The number of element in the list corresponds the length of
\code{lambda1}.}

\item{\code{lambda1}}{Vector (class \code{"numeric"}) of
\eqn{\ell_1}{l1} norm penalty levels for which the model has
eventually been fitted.}

\item{\code{lambda2}}{Scalar for the level of (structuring)
\eqn{\ell_2}{l2} norm penalty applied.}

\item{\code{df}}{Estimated degree of freedoms for the successive
\code{lambda1}.}

\item{\code{residuals}}{List of Matrices of residuals. Each entry
of the list corresponds to a value of \code{lambda1}.}
\item{\code{r.squared}:}{Vector (class \code{"numeric"}) given the
coefficient of determination as a function of lambda1.}

\item{\code{normx}}{Vector (class \code{"numeric"}) containing the
square root of the sum of squares of each column of the design
matrix.}

\item{\code{monitoring}}{List (class \code{"list"}) which
contains evolution of the loglikelihood and of the objective
function for each lambda1 and each iteration of the algorithm.}

\item{\code{fitted}}{List of Matrices of fitted values. Each
entry of the list corresponds to a value of \code{lambda1}.}
}}
\section{Methods}{
  This class comes with the usual \code{predict(object,
  newx, ...)}, \code{fitted(object, ...)},
  \code{residuals(object, ...)}, \code{print(object, ...)},
  \code{show(object)} and \code{deviance(object, ...)}
  generic (undocumented) methods.

  The two methods \code{loglik(object)} and
  \code{objective(object)} are also defined to recover the
  log-likelihood and objective functions along the fit.
  They can be applied on a spring object or on a list of
  spring object. The corresponding output is an object with
  class \code{criterion}. See the corrdesponding
  documentation for defintiion and methods.

  The method \code{spr.model(object, lambda1,
  lambda2=NULL))} applies on a spring object or a list of
  spring objects to extract an single object with class
  \code{model}.  See the corresponding documentation for
  definition and methods.
}
\seealso{
See also \code{\link{spring}},
\code{\linkS4class{criterion}}, \code{\linkS4class{model}}.
}
\keyword{class}

