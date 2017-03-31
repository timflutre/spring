rm(list=ls()) ; gc()
library(spring)
source("plot.R")

cat("\nPRETREATMENTS")
cat("\n==========================================================")

cat("\n- Load the Brassica Napus dataset (pretreatment done)...")
## contains markers, traits and the genetic map
load("bnapus.RData")

## Problem size
n <- nrow(traits)
p <- ncol(markers)
q <- ncol(traits)

dist <- "genetic" ##
dist <- "neighborhood"

if (dist == "genetic") {
  ## Build Linkage Desequilibrium using the genetic map
  cat("\n- Prior construction for structure inference (genetic-distance based)...")
  d <- as.numeric(unlist(with(map, tapply(loci, chrom, function(x) c(diff(x), Inf)) ))) ; d <- d[-p]
  U <- bandSparse(p,p,k=0:1,diagonals=list(rep(1,p), -exp(-2*d)))
  V <- Diagonal(x=c(1/(1-exp(-2*2*d)),1))
  L <- t(U) %*% V %*% U
} else {
  ## Neighborhood base... (nice one)
  cat("\n- Prior construction for structure inference (neighborhood based)...")
  D1 <- diffOp(p,4)
  L <- t(D1) %*% D1
}

L <- Diagonal(x=1/sqrt(diag(L))) %*% L %*% Diagonal(x=1/sqrt(diag(L)))
normx <- colSums(markers^2)/(n-1)
L <- Diagonal(x=normx) %*% L %*% Diagonal(x=normx)

cat("\n\nNOW TO THE MAIN FITTING")
cat("\n==========================================================")

cat("\nFix lambda2 (prior) on a grid of lambda1 (sparsity) \n")
out  <- spring(markers, traits, nlambda1=20, min.ratio=3e-1, struct=L, lambda2=c(10^seq(3,-3,len=7),0))

cat("\nMODEL CHOICE AND OUTPUT")
cat("\n==========================================================")
cat("\n- Compute BIC/AIC information criteria \n")
crit <- pen.crit(out)

spring.reg <- crit$BIC$model@coef.regres
spring.dir <- crit$BIC$model@coef.direct
cov.res    <- crit$BIC$model@covariance

map.reg <- plot.map( spring.reg, map, title="", plot=FALSE, show.names=TRUE, chr=c(2,8,10))
map.dir <- plot.map(-spring.dir, map, title="", plot=FALSE, show.names=TRUE, chr=c(2,8,10))

readline("hit a key for plotting AIC/BIC used for model selection...")
plot(crit$BIC$criterion)

readline("hit a key for plotting the coefficents...")

grid.arrange(
    arrangeGrob(
        map.reg + theme(legend.position="none",plot.margin= unit(c(0,0.05,0,0), "lines")),
        map.dir + theme(legend.position="none",plot.margin= unit(c(0,0.05,0,0), "lines")),
        g_legend(map.dir), nrow=1, widths=c(5,5, 1), main="")
    )

readline("hit a key for plotting the estimated covariance (renormalized to a correlation matrix)...")
plot(crit$BIC$model, type="covariance", corr=TRUE, legend=FALSE)
