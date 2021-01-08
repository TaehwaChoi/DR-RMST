###########################
# Last update : 08.17.2019    <= Please check!
###########################
# L = 10, 20
# Censor rate = 30%
###########################

rm(list=ls())
setwd("/Users/flowerlee0428/Dropbox/DR_RMST(2019)/R/Simulations/res_ver.0815")
getwd()

####################################
#### Generate simulation dataset####
####################################
simdata = function(n, tau) 
{
  # Covariates (X's)
    x1 = runif(n, -1, 1)
    x2 = runif(n, -1, 1)
    x3 = runif(n, -1, 1)
    x4 = runif(n, -1, 1)
    x5 = runif(n, -1, 1)
  # expit 
    expit = function(x) (1 / (1+exp(-x)))
  # trt indicator -- X1 - X3 are involved
    A = rbinom(n, 1, prob = expit(x1 + x2*x3)) 
  # True parameters 
    model = exp(-5 - abs(x1) + 10*A*x2*x3)     # rate model of true surv time  
    Time  = rexp(n, rate = model)              # (unknown) true time 
  # Given information (our dataset)   
    Cens  = rexp(n, tau)                       # Censoring rate ~ 30%         
    Y     = pmin(Time, Cens)                   # Observed time (right censored)
    delta = ifelse(Time <= Cens, 1, 0)         # mean(I(T<=C)) = Prob(T<=C)            
    range(Y[A==0]) ; range(Y[A==1])
    1-mean(delta)
    dat = data.frame(Y = round(Y,5), delta, A, x1, x2, x3, x4, x5)
}

####################################
#### Results                    ####
####################################
Result = function(parest, sdest, truepar)   
{
  True = rep(truepar, ncol(parest))
  Est   = colMeans(parest)
  Bias  = Est - True 
  SEE   = apply(parest, 2, sd)
  ASE   = apply(sdest,  2, mean)
  ci.lower = t(parest) - (1.96 * t(sdest)) < True
  ci.upper = t(parest) + (1.96 * t(sdest)) > True
  Cover    = apply(ci.lower*ci.upper, 1, mean)
  result = rbind(Est, Bias, ESD = SEE, ASE = ASE, Cover)
}








