##' Class "model"
##'
##' Class of object returned by the \code{criteria} or
##' \code{cv.spring} functions.
##'
##' @aliases plot,model-method
##'
##' @docType class
##'
##' @keywords class
##'
##' @name model-class
##' @rdname model-class
##'
##' @exportClass model
##' @exportMethod plot
##'
##' @importFrom stats predict residuals deviance
##' @import ggplot2
##' @import reshape2
##' @import gridExtra
##' 
setClassUnion("pen", c("NULL","numeric"))
setClass("model",
  slots = list(
  lambda1     = "pen"   ,
  lambda2     = "pen"   ,
  coef.direct = "Matrix"    ,
  coef.regres = "Matrix"    ,
  covariance  = "Matrix"    ,
  intercept   = "Matrix")
)

setMethod("predict", "model", definition =
   function (object, newx, ...)  {
     if (is.null(object@intercept)) {
       return(newx %*% object@coef.regres)
     } else {
       return(sweep(newx %*% object@coef.regres,2L,-object@intercept,check.margin=FALSE))
     }
   }
)

setMethod("residuals", "model", definition =
   function (object, newx, newy, ...)  {
     return(newy - predict(object, newx))
   }
)

setMethod("deviance", "model", definition =
   function (object, newx, newy, ...)  {
     return(colSums(residuals(object, newx, newy)^2))
   }
)

setMethod("plot", "model", definition =
    function(x, y, type=c("coefficients", "covariance"), coef.shape=c("matrix", "line"), pred.subset=1:p, resp.subset=1:q, plot=TRUE, corr=FALSE, legend=TRUE, same.scale=FALSE, ...) {

      p <- nrow(x@coef.direct)
      q <- ncol(x@coef.direct)
      
      type <- match.arg(type)
      coef.shape <- match.arg(coef.shape)
      simple.plot <- switch(coef.shape, "matrix" = simple.matrix.plot, "line" = simple.line.plot)

      if (type %in% c("covariance")) {
        if (corr == TRUE) {
          cov <- cov2cor(as.matrix(x@covariance))
        } else {
          cov <- as.matrix(x@covariance)
        }
        d <- simple.matrix.plot(melt(cov)) +
          labs(x="outcome", y="outcome"  , title="Residual covariance")
        if (!legend) {
          d <- d + theme(legend.position="none")
        }
      } else {
        coef.direct <- x@coef.direct[pred.subset, resp.subset]
        coef.regres <- x@coef.regres[pred.subset, resp.subset]
        if (same.scale) {
          vmax <- max( max(coef.direct), max(coef.regres) )
          vmin <- min( min(coef.direct), min(coef.regres) )
          d.regres <- simple.plot(melt(as.matrix(coef.regres)), vmin, vmax) +
            labs(x="outcome", y="predictor", title="Regression coefficients")
          d.direct <- simple.plot(melt(as.matrix(coef.direct)), vmin, vmax) +
            labs(x="outcome", y="predictor", title="Direct effects")
        } else {
          d.regres <- simple.plot(melt(as.matrix(coef.regres))) +
            labs(x="outcome", y="predictor", title="Regression coefficients")
          d.direct <- simple.plot(melt(as.matrix(coef.direct))) +
            labs(x="outcome", y="predictor", title="Direct effects")
        }
        if (legend) {
          d.regres <- d.regres + theme(plot.margin= unit(c(0,0.05,0,0), "lines"))
          d.direct <- d.direct + theme(plot.margin= unit(c(0,0.05,0,0), "lines"))          
        } else {
          d.regres <- d.regres + theme(legend.position="none",plot.margin= unit(c(0,0.05,0,0), "lines"))
          d.direct <- d.direct + theme(legend.position="none",plot.margin= unit(c(0,0.05,0,0), "lines"))          
        }

        d <- grid.arrange(arrangeGrob(d.regres, d.direct, nrow=1),main="")

        ## if (coef.shape == "matrix") {
        ##   d <- grid.arrange(arrangeGrob(d.regres, d.direct, nrow=1), g_legend(d.regres),main="")
        ## } else {
        ##   d <- grid.arrange(
        ##          arrangeGrob(
        ##            d.regres + theme(legend.position="none",plot.margin= unit(c(0,0.05,0,0), "lines")) +
        ##            labs(y="outcome", x="predictor", title="Regression coefficients"),
        ##            d.direct + theme(legend.position="none",plot.margin= unit(c(0,0.05,0,0), "lines")) +
        ##            labs(y="outcome", x="predictor", title="Direct effects"),
        ##            nrow=2), g_legend(d.regres),  nrow=1, widths=c(8,2), main="")
        ## }
      }
      
      
      if (plot) {print(d)}

      return(invisible(d))
      
    })

simple.line.plot <- function(dplot, vmin=-max(abs(dplot$value)), vmax=max(abs(dplot$value))) {
  
  return(ggplot(data=dplot,aes(Var1,value,colour=as.factor(Var2),group=as.factor(Var2))) +
         geom_line() + geom_point() + ylim(c(vmin,vmax)) + scale_colour_discrete(name = "Outcome"))
  
}

simple.matrix.plot <- function(dplot, vmin=-max(abs(dplot$value)), vmax=max(abs(dplot$value))) {

  return(ggplot(dplot, aes(Var2, Var1, fill = value)) + geom_tile(colour="gray80") +
         scale_fill_gradient2(name="value",limits=c(vmin,vmax),low = muted("blue"),
                              mid = "white", high = muted("red"),
                              midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar") +
         coord_fixed() + theme_bw() + theme(legend.title  = element_text(size = 12),
                                            axis.text  = element_blank(),
                                            plot.background = element_blank(),
                                            panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(),
                                            panel.border = element_blank(),
                                            panel.background = element_blank()))
}

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

