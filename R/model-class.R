#' A Reference Class to represent a model in spring
#'
#' @description Class of object returned by the \code{criteria} or \code{cv.spring} functions.
#'
#' @field 
#'  
#' @import ggplot2
#' @import Matrix
#' @import reshape2
#' @import grid
#' @import gridExtra
#' @exportClass model
#' @export model
model <- setRefClass("model",
  field = list(
  coef.direct = "Matrix",
  coef.regres = "Matrix",
  covariance  = "Matrix",
  intercept   = "Matrix",
  lambda1 = "numeric",
  lambda2 = "numeric")
)

model$methods(predict = function(newx)  {
   return(sweep(newx %*% coef.regres, 2L, intercept, check.margin=FALSE))
})

model$methods(residuals = function (newx, newy)  {
     return(newy - .self$predict(newx))
})

model$methods(deviance = function (newx, newy)  {
     return(colSums(.self$residuals(newx, newy)^2))
})

model$methods(drawResponse = function(design) {
  
  stopifnot(ncol(design) == nrow(.self$coef.regres))
  
  n <- nrow(design)
  q <- ncol(.self$covariance)
  
  ## The stochastic noise
  noise <- as.matrix(matrix(rnorm(q*n),n,q) %*% chol(.self$covariance))
  
  ## The multivariate response matrix with covariance,
  ## heteroscedasticity and independent
  return(.self$predict(design) + noise)
})

model$methods(plot = function(type=c("coefficients", "covariance"), coef.shape=c("matrix", "line"), pred.subset=1:p, resp.subset=1:q, plot=TRUE, corr=FALSE, legend=TRUE, same.scale=FALSE) {

      p <- nrow(.self$coef.direct)
      q <- ncol(.self$coef.direct)
      
      type <- match.arg(type)
      coef.shape <- match.arg(coef.shape)
      simple.plot <- switch(coef.shape, "matrix" = simple.matrix.plot, "line" = simple.line.plot)

      if (type %in% c("covariance")) {
        if (corr == TRUE) {
          cov <- cov2cor(as.matrix(.self$covariance))
        } else {
          cov <- as.matrix(.self$covariance)
        }
        d <- simple.matrix.plot(melt(cov)) +
          labs(x="outcome", y="outcome"  , title="Residual covariance")
        if (!legend) {
          d <- d + theme(legend.position="none")
        }
      } else {
        p.coef.direct <- .self$coef.direct[pred.subset, resp.subset]
        p.coef.regres <- .self$coef.regres[pred.subset, resp.subset]
        if (same.scale) {
          vmax <- max( max(p.coef.direct), max(p.coef.regres) )
          vmin <- min( min(p.coef.direct), min(p.coef.regres) )
          d.regres <- simple.plot(melt(as.matrix(p.coef.regres)), vmin, vmax) +
            labs(x="outcome", y="predictor", title="Regression coefficients")
          d.direct <- simple.plot(melt(as.matrix(p.coef.direct)), vmin, vmax) +
            labs(x="outcome", y="predictor", title="Direct effects")
        } else {
          d.regres <- simple.plot(melt(as.matrix(p.coef.regres))) +
            labs(x="outcome", y="predictor", title="Regression coefficients")
          d.direct <- simple.plot(melt(as.matrix(p.coef.direct))) +
            labs(x="outcome", y="predictor", title="Direct effects")
        }
        if (legend) {
          d.regres <- d.regres + theme(plot.margin= unit(c(0,0.05,0,0), "lines"))
          d.direct <- d.direct + theme(plot.margin= unit(c(0,0.05,0,0), "lines"))          
        } else {
          d.regres <- d.regres + theme(legend.position="none",plot.margin= unit(c(0,0.05,0,0), "lines"))
          d.direct <- d.direct + theme(legend.position="none",plot.margin= unit(c(0,0.05,0,0), "lines"))          
        }

        d <- grid_arrange_shared_legend(d.regres, d.direct, nrow=1, ncol=2)
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

grid_arrange_shared_legend <- function(..., nrow = 1, ncol = length(list(...)), position = c("bottom", "right")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.draw(combined)
}
