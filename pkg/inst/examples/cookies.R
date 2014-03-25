rm(list=ls()); gc()
library(spring)
library(parallel)
source("plot.R")

cat("\nPRETREATMENTS")
cat("\n==========================================================")

cat("\n- Load the cookies dataset (pretreatment done)...")
load("cookies.RData")

## Problem dimension
n <- nrow(chem.comp)
p <- ncol(nri)
q <- ncol(chem.comp)

## let us try oth order derivate operator
o <- 1
cat("\n- Prior construction for structure inference (neighborhood based)...")
D <- diffOp(p,o)
L <- t(D) %*% D
L <- Diagonal(x=1/sqrt(diag(L))) %*% L %*% Diagonal(x=1/sqrt(diag(L)))
normx <- colSums(nri^2)/(n-1)
L <- Diagonal(x=normx) %*% L %*% Diagonal(x=normx)

train <- 1:39
test <- setdiff(1:n,train)

cat("\n\nNOW TO THE MAIN FITTING")
cat("\n==========================================================")
sp.ic <- spring(nri[train, ], chem.comp[train, ], nlambda1=20, min.ratio=5e-4, struct=L, lambda2=10^seq(-1,-4,len=4))

cat("\n\nCROSS-VALIDATION")
cat("\n==========================================================")
sp.cv <- cv.spring(nri[train, ], chem.comp[train, ], K=4, nlambda1=10, min.ratio=5e-4, struct=L, lambda2=10^seq(-1,-4,len=4))

cat("\nMODEL CHOICE AND OUTPUT")
cat("\n==========================================================")
cat("\n- Compute BIC/AIC information criteria and get CV-min model \n")

BIC <- pen.crit(sp.ic)$BIC
AIC <- pen.crit(sp.ic)$AIC
CV  <- pen.crit(sp.cv)$cv.min

err.BIC <- colMeans(residuals(BIC$model, as.matrix(nri[test, ]), as.matrix(chem.comp[test, ]))^2)
err.AIC <- colMeans(residuals(AIC$model, as.matrix(nri[test, ]), as.matrix(chem.comp[test, ]))^2)
err.CV  <- colMeans(residuals(CV$model , nri[test, ], chem.comp[test, ])^2)

readline("hit a key for plotting BIC used for model selection...")
plot(BIC$criterion)

readline("hit a key for plotting the coefficents...")
plot(BIC$model, coef.shape="line")

readline("hit a key for plotting the estimated covariance...")
cov.plot(BIC$model@covariance, corr=FALSE)
