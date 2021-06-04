## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("factormodel")

## ----gh-installation, eval = FALSE--------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("yujunghwang/factormodel")

## -----------------------------------------------------------------------------
library(factormodel)
library(nnet)
library(pracma)
library(stats)
library(utils)

# DGP
# set parameters
nsam <- 5000

M1 <- rbind(c(0.8,0.1,0.1),c(0.1,0.2,0.7))
M2 <- rbind(c(0.7,0.2,0.1),c(0.2,0.2,0.6))
M3 <- rbind(c(0.9,0.05,0.05),c(0.1,0.1,0.8))

CM1 <- t(apply(M1,1,cumsum))
CM2 <- t(apply(M2,1,cumsum))
CM3 <- t(apply(M3,1,cumsum))

# 40% of sample is type 1, 60% is type 2
truetype <- as.integer(runif(nsam)<=0.4) +1

# generate fake data
dat <- data.frame(msr1=rep(NA,nsam),msr2=rep(NA,nsam),msr3=rep(NA,nsam))

for (k in 1:nsam){
  dat$msr1[k] <- which(runif(1)<=CM1[truetype[k],])[1]
  dat$msr2[k] <- which(runif(1)<=CM2[truetype[k],])[1]
  dat$msr3[k] <- which(runif(1)<=CM3[truetype[k],])[1]
}

# estimate using dproxyme
oout <- dproxyme(dat=dat,sbar=2,initvar=1,initvec=NULL,seed=210313,tol=0.005,maxiter=200,miniter=10,minobs=100,maxiter2=1000,trace=FALSE)

# check whether the estimated measurement stochastic matrices are same with the true # measurement stochastic matrices
print(oout$M_param)

# check type probability
print(head(oout$typeprob))

## -----------------------------------------------------------------------------
library(factormodel)
library(stats)
library(utils)
library(gtools)

set.seed(seed=210315)

# DGP
# set parameters
nsam <- 5000 # number of observations
np <- 3 # number of proxies

true_mtheta <- 2
true_vartheta <- 1.5
true_theta <- rnorm(nsam, mean=true_mtheta, sd=sqrt(true_vartheta))

# first proxy variable is an anchoring variable
true_alpha0 <- c(0,2,5)
true_alpha1 <- c(1,0.5,2)
true_varnu  <- c(0.5,2,1)

# simulate fake data
dat <- matrix(NA,nrow=nsam,ncol=np)
for (k in 1:np){
  dat[,k] <- true_alpha0[k] + true_alpha1[k]*true_theta + rnorm(nsam,mean=0,sd=sqrt(true_varnu[k]))
}

# estimate parameters using cproxyme
oout <- cproxyme(dat=dat,anchor=1)

# print estimated parameters
print(oout$alpha0)
print(oout$alpha1)
print(oout$varnu)
print(oout$mtheta)
print(oout$vartheta)

