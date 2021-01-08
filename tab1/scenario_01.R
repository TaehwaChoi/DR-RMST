##################################
## Scenario 1
## Model:: PS=Correct; OR=Correct
##################################

# rm(list=ls())

#### simulation ####
L <- 20       #tmax (setting | L=10; L=20)
n <- 600      #sample size 600 (600 random samples per boot)
tau <- exp(-4.5)
nsim <- 2000    #nsim=2000
nboot <- 200    # boot 200
truepar <- ifelse(L==10, 0.871, 1.682) 
Amod = "A ~ Z1 + Z2 + Z3"
Qmod = "prmst ~ Z1 + Z2 + Z3"

parest <- sdest <- NULL

for(i in 1:nsim) 
{
  set.seed(i)
  dt <- simdata(n, tau) 
  Y <- dt$Y
  delta <- dt$delta
  
  # pseudo-RMST (prmst)
  prmst <- dt$prmst <- pseudomean(time=Y, event=delta, tmax=L) #prmst: pseudo rmst
  
  # logistic reg model for IPW (TRUE)
  # mod.p stands for PS model 
  mod.p <- glm(Amod, data = dt, family = binomial)
  p.1 <- predict(mod.p, type = "response")  # p.1 = Prob(A=1|Z)
  p.0 <- (1-p.1)                            # p.0 = Prob(A=0|Z)

  #### Pseudo-obs-based IPW (Re-weigthed sample)
  df.est.ipw <- matrix(0,nrow=n, ncol=2)
  df.est.ipw[,1] <- (dt$prmst*(1-dt$A))/p.0   #est.ipw.0
  df.est.ipw[,2] <- (dt$prmst*dt$A)/p.1       #est.ipw.1
  m.ipw <- colMeans(df.est.ipw)
  
  # OR model (i.e., Q model) (CORRECT) 
  #### model-based direct adjustment (g computation) (true models)  
  mod.m0 <- lm(Qmod, data = subset(dt,A==0))
  mod.m1 <- lm(Qmod, data = subset(dt,A==1))
  
  df.est.m <- matrix(0,nrow=n, ncol=2)
  df.est.m[,1] <- predict(mod.m0, newdata = dt)
  df.est.m[,2] <- predict(mod.m1, newdata = dt) 
  m.qmod <- colMeans(df.est.m)
  
  #### DR: yi = Di*yi(1) + (1-Di)*yi(0) 
  df.est.dr <- matrix(0,nrow=n, ncol=2)
  # DR estimator for E(yi(0))
  df.est.dr[,1] <- (prmst*(1-dt$A))/p.0 - (((1-dt$A) - p.0)*df.est.m[,1])/p.0
  # DR estimator for E(yi(1))
  df.est.dr[,2] <- (prmst*dt$A)/p.1 - ((dt$A - p.1)*df.est.m[,2])/p.1
  m.dr <- colMeans(df.est.dr)
  
  # result table of ith simulation 
  est <- matrix(c(m.ipw,0,m.qmod,0,m.dr,0), nrow=3, ncol=3, byrow=T)
  rownames(est) <- c("ipw","qmod","dr"); colnames(est) <- c("A==0","A==1","diff")
  est[,3]<-c(est[,2]-est[,1])
  parest <- rbind(parest,c(est[,3]))
  
  #### bootstrap ####
  parest.boot <- NULL 
  
  for(b in 1:nboot){ #nboot=60
    
    idx <- sample(1:n,n,replace=T)
    dt.boot <- dt[idx,] #dt should not be changed
    
    # pseudo-RMST (prmst)
    prmst.boot <- dt.boot$prmst <- pseudomean(time=dt.boot$Y, event=dt.boot$delta, tmax=L) #prmst: pseudo rmst
    
    # logistic reg model for IPW (TRUE)
    # mod.p stands for PS model 
    mod.p <- glm(Amod, data = dt.boot, family = binomial)
    p.1.boot <- predict(mod.p, type = "response")       # p.1 = Prob(A=1|Z)
    p.0.boot <- (1-p.1.boot)                            # p.0 = Prob(A=0|Z)
    
    #### Pseudo-obs-based IPW (Re-weigthed sample)
    df.ipw.boot <- matrix(0,nrow=n, ncol=2)
    df.ipw.boot[,1] <- (prmst.boot)*((1-dt.boot$A)/p.0.boot) #est.ipw.0
    df.ipw.boot[,2] <- (prmst.boot)*(dt.boot$A/p.1.boot)     #est.ipw.1
    m.ipw.boot <- colMeans(df.ipw.boot)
    
    # OR model (i.e., Q model) (TRUE) 
    #### model-based direct adjustment (true models)  
    mod.m0 <- lm(Qmod, data = subset(dt.boot,A==0))
    mod.m1 <- lm(Qmod, data = subset(dt.boot,A==1))
    
    df.q.boot <- matrix(0,nrow=n, ncol=2)
    df.q.boot[,1] <- predict(mod.m0, newdata = dt.boot)
    df.q.boot[,2] <- predict(mod.m1, newdata = dt.boot) 
    m.q.boot <- colMeans(df.q.boot)

    #### DR: yi = Di*yi(1) + (1-Di)*yi(0) 
    df.dr.boot <- matrix(0,nrow=n, ncol=2)
    # DR estimator for E(yi(0))
    df.dr.boot[,1] <- ((prmst.boot * (1-dt.boot$A) ) / p.0.boot) - ((( (1-dt.boot$A) - p.0.boot)*df.q.boot[,1] ) / p.0.boot)
    # DR estimator for E(yi(1))
    df.dr.boot[,2] <- (( prmst.boot * dt.boot$A ) / p.1.boot) - (( (dt.boot$A - p.1.boot)*df.q.boot[,2] ) / p.1.boot)
    m.dr.boot <- colMeans(df.dr.boot)

    # result table of ith simulation 
    est.boot <- matrix(c(m.ipw.boot,0,m.q.boot,0,m.dr.boot,0), nrow=3, ncol=3, byrow=T)
    rownames(est.boot) <- c("ipw","qmod","dr"); colnames(est.boot) <- c("A==0","A==1","diff")
    
    est.boot[,3]<-c(est.boot[,2]-est.boot[,1])
    parest.boot <- rbind(parest.boot,c(est.boot[,3]))
  } # Boot-loop ends here. 

   sdest <- rbind(sdest, apply(parest.boot,2,sd))

   cat("========== i =",i,": Scenario 1 ===========","\n")

     if(i>1){
     result = Result(parest,sdest,truepar)
     print(result)

     }
   
   flname=paste("result_",n,"_",L,"_scenario1.txt",sep="")
   if(i%%10==0) write.table(result,file=flname)
}
 



