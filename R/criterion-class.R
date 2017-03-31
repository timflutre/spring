##' Class "criterion"
##'
##' Class of object returned by the \code{criteria} or
##' \code{cv.spring} functions.
##'
##' @slot name a character describing the criterion (AIC, BIC,
##' loglikelihood, objective function, CV.min, CV.1se)
##'
##' @slot data a data frame with evolution of the criterion as a
##' function of \code{lambda1} and \code{lambda2} and the
##' corresponding degrees of freedom and shrinkage factor.
##'
##' @docType class
##'
##' @keywords class
##'
##' @seealso See also \code{\link{plot,criterion-method}}.
##'
##' @name criterion-class
##' @rdname criterion-class
##'
##' @exportClass criterion
##'
setClass("criterion",
  slots = list(
  name        = "character" ,
  data        = "data.frame")
)

##' Plot method for object criterion
##'
##' A plotting method for the class \code{criterion}.
##'
##' @usage \\S4method{plot}{criterion}(x, y, xvar="lambda1",
##' log.scale=ifelse(xvar=="lambda1",TRUE,FALSE), plot=TRUE, ...)
##'
##' @param x an object with class criterion
##' @param y used for S4 compatibility.
##' @param xvar a character among \code{"lambda1"}, \code{"fraction"}
##' or \code{"df"} for the x-axis (only \code{lambda1} is available
##' when the criterion object come from a call to \code{cv.spring}).
##' @param log.scale a boolean to put the x-axis on a log-scale.
##' @param plot is \code{TRUE}, a plot is produce. Otherwise a
##' \code{ggplot} object is sent back.
##'
##' @return a invisible ggplot object.
##'
##' @seealso \code{\linkS4class{criterion}}.
##'
##' @name plot,criterion-method
##' @aliases plot,criterion-method plot.criterion
##' @docType methods
##' @rdname plot.criterion
##'
##' @import ggplot2 scales grid reshape2
##' @exportMethod plot
##' @export
setMethod("plot", "criterion", definition =
  function(x, y, xvar="lambda1", log.scale=ifelse(xvar=="lambda1",TRUE,FALSE), plot=TRUE, ...) {

    x@data$lambda2 <- as.factor(x@data$lambda2)
    x@data$xvar <- switch(xvar,
                          "lambda1"  = x@data$lambda1,
                          "df"       = x@data$df     ,
                          "fraction" = x@data$fraction)
    d <- ggplot(x@data, aes(x=xvar, y=value, colour=lambda2, group=lambda2))

    xlab <- switch(xvar,
                   "fraction" = expression(paste("|",beta[lambda[1]],"|",{}[1]/max[lambda[1]],"|",beta[lambda[1]],"|",{}[1],sep="")),
                   "df" = "Estimated degrees of freedom",
                   ifelse(log.scale,expression(log[10](lambda[1])),expression(lambda[1])) )

    d <- d +  geom_line(aes(x=xvar,y=value)) + geom_point(aes(x=xvar,y=value)) +
        labs(x=xlab, y="criterion's value",  title=x@name)

    ## if (!is.null(x@data$serr)) {
    ##   d <- d + geom_smooth(aes(ymin=value-serr, ymax=value+serr), data=x@data, alpha=0.2, stat="identity")
    ## }

    d <- d + scale_colour_discrete(name = expression(lambda[2]))
    d <- d + theme(legend.title = element_text(size = 12))

    if (log.scale) {
      d <- d + scale_x_log10() + annotation_logticks(sides="b")
    }

    if (plot) {print(d)}

    return(invisible(d))

  })
