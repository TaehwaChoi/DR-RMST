#### Simulation Setting ###############################################
### 1. Baseline covariate vector Z = {Z1, Z2, Z3}' 
###    : multivariate normal with mean zero, unit var, corr(Z1,Z3)=0.2
### 2. trt indicator A ~ Bernoulli; p=expit(-0.5Z1-0.5*Z2))
### 3. T from an exponential dist. 
###    : for A=0, par = exp(-2.5-1.5*Z1-1*Z2-0.7*Z3)
###    : for A=1, par = exp(-3-1*Z1-0.9*Z2-1*Z3) 
### 4. censoring time C ~ exponential with exp(-4.5)
####   which lead to app. 25 % censoring 
#######################################################################
### last update : 05 April  

library(pseudo); library(survival); library(KMsurv)
library(mvtnorm);library(geepack); library(ggplot2)
getwd() # "/Users/flowerlee0428/Desktop"

  
#####################################
#### Generate simulation dataset####
####################################

rm(list=ls())
### data 
simdata <- function(n, tau) { 
  Sigma <- diag(rep(1,3)); Sigma[1,3]=Sigma[3,1]=0.2 
  Z <- rmvnorm(n, mean=rep(0,3), sigma=Sigma, method="chol")
  Z1 <- Z[,1]; Z2 <- Z[,2]; Z3 <- Z[,3] 
  # trt indicator (Bernoulli)
  A <- rbinom(n,1,expit(-0.5*Z1-0.5*Z2))# Z1 and Z2 are involved in trt assignment 
  par.t <- exp((-2.5-1.5*Z1-1*Z2-0.7*Z3)*(A==0) + (-3-1*Z1-0.9*Z2-1*Z3)*(A==1)) 
  T <- rexp(n, par.t)                   # true surv time 
  C <- rexp(n, tau)                     # independent censoring 
                                          # tau as a controller of censoring rate 
  Y <- pmin(T,C)                        # observed time 
  delta <- ifelse(T<=C,1,0)             # censoring indicator  
  simdata <- data.frame(Y = round(Y,5),delta,A,Z1,Z2,Z3,T)
  simdata[order(Y),]
}

### function 
expit <- function(x) exp(x)/(1+exp(x))

#### result #### **CHECK AGAIN! 
Result <- function(parest, sdest, turepar) { #truepar <- L 
  param <- rep(truepar, ncol(parest))
  Est <- colMeans(parest)
  Bias <- Est-param
  ASE <- apply(sdest,2,mean) #avg of 1500 estimates 
  SEE <- apply(parest,2,sd)  #se estimate 
  ci.lower = t(parest) - (1.96 * t(sdest)) < param 
  ci.upper = t(parest) + (1.96 * t(sdest)) > param  
  Cover <- apply(ci.lower*ci.upper, 1, mean)
  #res.tab <- rbind(Est, Bias, ASE, SEE, CI)
  result <- rbind(Est, Bias, ASE, SEE, Cover)
}
