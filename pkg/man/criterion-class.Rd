\docType{class}
\name{criterion-class}
\alias{criterion-class}
\title{Class "criterion"}
\description{
Class of object returned by the \code{criteria} or
\code{cv.spring} functions.
}
\section{Slots}{

\describe{
\item{\code{name}}{a character describing the criterion (AIC, BIC,
loglikelihood, objective function, CV.min, CV.1se)}

\item{\code{data}}{a data frame with evolution of the criterion as a
function of \code{lambda1} and \code{lambda2} and the
corresponding degrees of freedom and shrinkage factor.}
}}
\seealso{
See also \code{\link{plot,criterion-method}}.
}
\keyword{class}

