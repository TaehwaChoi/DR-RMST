###############################
# Last update : 06.22.2019       <= Please check!
###############################
# Case1 : Q TRUE, A TRUE 
###############################
# TRUE (L=10) 
# RMST_A1   RMST_A0     delta 
# 8.430547  9.790229 -1.359682 
###############################
# L = 10, 20
# Censor rate = 30%
###############################
rm(list=ls())

#### Load functions : simdata, Result 
source("/Users/flowerlee0428/Dropbox/DR_RMST(2019)/R/Simulations/res_ver.0815/sim_dat0815.R")
#### simulation ##############################################
CASE  = 1
L     = 10          # tmax 
n     = 600         # sample size 
nsim  = 2000        # number of iteration 
#nsim  = 100
nboot = 200          # sample size per bootstrap (of each iter)
tau  = exp(-6.7)    # determines censoring rate  
truepar = ifelse(L==10, -1.360, -3.186)
###############################################################
#### Load packages 
{
  library(pseudo); library(survival); library(KMsurv)
  library(mvtnorm);library(geepack); library(ggplot2)
  library(ranger); library(survRM2); library(xtable); library(mgcv)
  library(progress)
}

#### Progress bar 
ProgressBar = progress_bar$new(total = nsim)

parest  =  NULL
sdest   =  NULL

for (i in 1:nsim) { 
  
  ProgressBar$tick()
  
  ################################ estimates ################################
  # Generate dataset
  set.seed(i)
  dt     =  simdata(n, tau)
  Y      =  dt$Y
  delta  =  dt$delta
  
  # Pseudo mean value generation
  dt$prmst = prmst = pseudomean(time = Y, event = delta, tmax = L)
  
  # Q mod (OLS) = TRUE 
  lmfit1 = lm(prmst ~ x1 + x2 + x3 + x4 + x5, data = subset(dt, A==1))
  lmfit0 = lm(prmst ~ x1 + x2 + x3 + x4 + x5, data = subset(dt, A==0))
  lm.Q1  = predict(lmfit1, newdata = dt)         # yi(1), i = 1, ..., n 
  lm.Q0  = predict(lmfit0, newdata = dt)         # yi(0), i = 1, ..., n
  
  # A mod (Logistic Reg) = TRUE  
  lm.A = glm(A ~ x1 + x2 + x3 + x4 + x5, family = "binomial", data = dt)
  pi1 = predict(lm.A, type = "response")         # pi1 = Prob(A=1|X)
  pi0 = (1 - pi1)                                # pi0 = Prob(A=0|X)
  lm.A1 = (dt$prmst *    dt$A)  / pi1            # Confounding effect correction 
  lm.A0 = (dt$prmst * (1-dt$A)) / pi0            # Confounding effect correction
  
  # Q mod (RF) = TRUE
  rffit1 = ranger(prmst ~ x1 + x2 + x3 + x4 + x5, data = dt[dt$A==1,] )
  rffit0 = ranger(prmst ~ x1 + x2 + x3 + x4 + x5, data = dt[dt$A==0,] )
  rf.Q1  = predict(rffit1, data = dt)$prediction
  rf.Q0  = predict(rffit0, data = dt)$prediction
  
  ## DR: yi = Di*yi(1) + (1-Di)*yi(0) 
  # DR (lin=TRUE, lin=TRUE) 
  # DR (lin, lin) 
  DR.m1 = (prmst *    dt$A)  /  pi1 - ((   dt$A  - pi1) * lm.Q1) / pi1          # DR estimate for yi(1)
  DR.m0 = (prmst * (1-dt$A)) /  pi0 - (((1-dt$A) - pi0) * lm.Q0) / pi0          # DR estimate for yi(0)
  
  # DR (RF = TRUE, logistic = TRUE)
  # DR (RF, Logistic) / col 16:18
  DRrf.m1 = (prmst *    dt$A)  /  pi1 -  ((   dt$A   -  pi1) * rf.Q1) / pi1  # DR estimate for yi(1)
  DRrf.m0 = (prmst * (1-dt$A)) /  pi0 -  (((1-dt$A)  -  pi0) * rf.Q0) / pi0  # DR estimate for yi(0)
  
  # diff (hat_delta)
  est_delta.Qlin  =  mean(lm.Q1) -  mean(lm.Q0) 
  est_delta.Alin  =  mean(lm.A1) -  mean(lm.A0)
  est_delta.Qrf   =  mean(rf.Q1) -  mean(rf.Q0)
  est_delta.DRlin =  mean(DR.m1) -  mean(DR.m0)
  est_delta.DRrf  =  mean(DRrf.m1) - mean(DRrf.m0)

  # result table of i-th iteration 
  parest = rbind(parest, c(est_delta.Qlin,  est_delta.Alin, 
                           est_delta.Qrf,   
                           est_delta.DRlin, est_delta.DRrf))
  colnames(parest) = c("Q(lin)","A(lin)","Q(RF)","DR(lin)","DR(RF)")
  
  ################################ bootstrap ################################
  
  parest_boot = NULL
  
  for(b in 1:nboot) 
    {
    idx = sample(1:n, n, replace = TRUE)
    dt.boot = dt[idx, ]
    
    # pseudo-rmst 
    prmst =  dt.boot$prmst = pseudomean(time=dt.boot$Y, event = dt.boot$delta, tmax = L)
     
    # Q mod (OLS) = TRUE 
    # Q mod (lin) / col 1:3x
    lmfit1 = lm(prmst ~ x1 + x2 + x3 + x4 + x5, data = subset(dt.boot, A==1))
    lmfit0 = lm(prmst ~ x1 + x2 + x3 + x4 + x5, data = subset(dt.boot, A==0))
    lm.Q1  = predict(lmfit1, newdata = dt.boot)         # yi(1), i = 1, ..., n 
    lm.Q0  = predict(lmfit0, newdata = dt.boot)         # yi(0), i = 1, ..., n
    
    # A mod (IPTW / Logistic Reg) = TRUE  
    # A mod (lin) / col 4:6
    lm.A = glm(A ~ x1 + x2 + x3 + x4 + x5, family = "binomial", data = dt.boot)
    pi1 = predict(lm.A, type = "response")                   # pi1 = Prob(A=1|X)
    pi0 = (1 - pi1)                                          # pi0 = Prob(A=0|X)
    lm.A1 = (dt.boot$prmst *    dt.boot$A)  / pi1            # Confounding effect correction 
    lm.A0 = (dt.boot$prmst * (1-dt.boot$A)) / pi0            # Confounding effect correction
    
    # Q mod (RF) = TRUE
    # Q mod (RF) / col 7:9
    rffit1 = ranger(prmst ~ x1 + x2 + x3 + x4 + x5, data = dt.boot[dt.boot$A==1,] )
    rffit0 = ranger(prmst ~ x1 + x2 + x3 + x4 + x5, data = dt.boot[dt.boot$A==0,] )
    rf.Q1  = predict(rffit1, data = dt.boot)$prediction
    rf.Q0  = predict(rffit0, data = dt.boot)$prediction
    
    ## DR: yi = Di*yi(1) + (1-Di)*yi(0) 
    # DR (lin=TRUE, lin=TRUE) 
    # DR (lin, lin) / col 13:15
    DR.m1 = (prmst *    dt.boot$A)  /  pi1 - ((   dt.boot$A  - pi1) * lm.Q1) / pi1       # DR estimate for yi(1)
    DR.m0 = (prmst * (1-dt.boot$A)) /  pi0 - (((1-dt.boot$A) - pi0) * lm.Q0) / pi0       # DR estimate for yi(0)
    
    # DR (RF = TRUE, RF=TRUE)
    # DR (RF, RF) / col 16:18
    DRrf.m1 = (prmst *    dt.boot$A)  /  pi1 -  ((   dt.boot$A   -  pi1) * rf.Q1) / pi1  # DR estimator for E(yi(1))
    DRrf.m0 = (prmst * (1-dt.boot$A)) /  pi0 -  (((1-dt.boot$A)  -  pi0) * rf.Q0) / pi0  # DR estimator for E(yi(0))
    
    # diff (hat_delta)
    est_delta.Qlin  =  mean(lm.Q1) -  mean(lm.Q0) 
    est_delta.Alin  =  mean(lm.A1) -  mean(lm.A0)
    est_delta.Qrf   =  mean(rf.Q1) -  mean(rf.Q0)
    est_delta.DRlin =  mean(DR.m1) -  mean(DR.m0)
    est_delta.DRrf  =  mean(DRrf.m1) - mean(DRrf.m0)
    
    
    # result table of i-th simulation's b-th bootstrap (rbind)
    parest_boot  = rbind(parest_boot, c(est_delta.Qlin,  est_delta.Alin, 
                                        est_delta.Qrf,  
                                        est_delta.DRlin, est_delta.DRrf ))
    colnames(parest_boot) = c("Q(lin)","A(lin)","Q(RF)","DR(lin)","DR(RF)")
    
    
    } # ----- bootstrap ends 
  
  ## Compute ESD by using bootstrap results  
  sdest = rbind(sdest, apply(parest_boot,2,sd))
    
  cat("========== i =",i,": Scenario 1 ===========","\n")
    
  if(i>1)
    {
      result = Result(round(parest,3), round(sdest,3), round(truepar,3))
      print(result)
    }
  
  flname = paste("Result_",n,"_",L,"_CASE",CASE,".txt",sep="")
  if(i%%10 == 0) write.table(result, file = flname)
  
  }

############################################################################################################

{
  save(result, file = "res_L10_CASE1.RData")
  save(sdest,  file = "res_sdest_L10_CASE1.RData")
  save(parest, file = "res_parest_L10_CASE1.RData")
}


