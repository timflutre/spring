##' Generate multivariate data 
##' 
##' @examples 
##' \dontrun{
##' ##  Basic simulation settings with no particular structure
##' k <- 10
##' p <- 20
##' q <- 3
##' ## direct coefficients
##' Omega.xy <- Matrix(0,p,q)
##' Omega.xy[sample(1:(p*q),k)] <- sample(c(-1,1),k,rep=T)
##' ## covariance of the noise (residual)
##' cor <- 0.5
##' R <- toeplitz(cor^(1:q-1))
##' ## turn this into a "spring" model
##' true.model <- new("model",
##'                   coef.direct = Omega.xy,
##'                   coef.regres = Matrix(- Omega.xy %*% R),
##'                   covariance  = Matrix(R),
##'                   intercept   = Matrix(rep(0,q),q,1))
##' 
##' ## draw some data with a given design matrix
##' dataset <- spring.generator(true.model, design)
##' }
##' 
##' @export
spring.generator <- function(model, design) {

  stopifnot(ncol(design) == nrow(model@coef.regres))

  n <- nrow(design)
  q <- ncol(model@covariance)
  
  ## The stochastic noise
  noise <- as.matrix(matrix(rnorm(q*n),n,q) %*% chol(model@covariance))
  
  ## The multivariate response matrix with covariance,
  ## heteroscedasticity and independent
  y <- predict(model, design) + noise
  r2 <- 1-colSums(noise^2)/colSums(y^2)
  
  return(list(design=design, response=y, model=model, r2=r2))
}

# x <- dataset$x
# y <- dataset$y
# train <- dataset$train
# test  <- dataset$test
