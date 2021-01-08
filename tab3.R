pkg = c("survival", "pseudo", "tmle", "MASS","mvtnorm")
sapply(pkg, library, character.only=T)
# test1: L = 10, TT
# test2: L = 10, TF (Z1)
# test3: L = 10, FT (Z1, Z3)
# test4: L = 10, FF

# test5: L = 20, TT
# test6: L = 20, TF
# test7: L = 20, FT
# test8: L = 20, FF
truepar = c(0.871, 1.682) 
expit = function(x) exp(x) / (1 + exp(x))
simdata1 <- function(n, tau) { 
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
n = 600; tau = exp(-4.5); L = 10
nsim = 2000
modelspec = function(Q = c("T","F"), A = c("T","F")){
   if(Q == "T" & A == "T") {
      Qmod = "prmst ~ Z1 + Z2 + Z3"
      Qmodt = "Y ~ A*(Z1 + Z2 + Z3)"
      Amod = "as.factor(A) ~ Z1 + Z2 + Z3"
   } else if(Q == "T" & A == "F") {
      Qmod = "prmst ~ Z1 + Z2 + Z3"
      Qmodt = "Y ~ A*(Z1 + Z2 + Z3)"
      Amod = "as.factor(A) ~ Z1"
   }  else if(Q == "F" & A == "T") {
      Qmod = "prmst ~ Z1 + Z3"
      Qmodt = "Y ~ A*(Z1 + Z3)"
      Amod = "as.factor(A) ~ Z1 + Z2 + Z3"
   } else {
      Qmod = "prmst ~ Z1 + Z3"
      Qmodt = "Y ~ A*(Z1 + Z3)"
      Amod = "as.factor(A) ~ Z1"
   }
   list(Qmod = Qmod, Qmodt = Qmodt, Amod = Amod)
}



test1 = NULL
ms = modelspec(Q = "T", A = "T")
Qmod = ms$Qmod; Qmodt = ms$Qmodt; Amod = ms$Amod


for (i in 1:nsim)
{
   set.seed(i)
   dt = simdata1(n, tau)
   Y = dt$Y; delta = dt$delta; A = dt$A; 
   Z1 = dt$Z1; Z2 = dt$Z2; Z3 = dt$Z3
   Xt = with(dt, as.matrix(data.frame(Z1, Z2, Z3)))
   
   
   dt$prmst = prmst = pseudomean(time = Y, event = delta, tmax = L)
   # dr
   m1 = predict(lm(Qmod, data = subset(dt, A == 1)), newdata = dt)
   m0 = predict(lm(Qmod, data = subset(dt, A == 0)), newdata = dt)
   pi1 = fitted(glm(Amod, data = dt, family = binomial))
   dr = mean((A/pi1 - (1-A)/(1-pi1))*prmst - 
                (A-pi1)/pi1*m1 - (A-pi1)/(1-pi1)*m0)
   
   
   # dr+
   idx1 = which(A==1); idx0 = which(A==0)
   d = log(mean(A*(1+exp(pi1)) / exp(pi1)))
   ENpi1 = as.numeric(exp(d + pi1) / (1+exp(pi1)))
   ENpi1 = pi1
   wt1 = ((1-ENpi1)/ENpi1^2)[idx1]
   wt0 = ((ENpi1)/(1-ENpi1)^2)[idx0]
   EN1 = predict(lm(Qmod, data = subset(dt, A ==1), weights = wt1), dt)
   EN0 = predict(lm(Qmod, data = subset(dt, A ==0), weights = wt0), dt)
   endr = mean((A/ENpi1-(1-A)/(1-ENpi1))*prmst - (A-ENpi1)/ENpi1*EN1 - (A-ENpi1)/(1-ENpi1)*EN0)
   
   
   # TMLE
   # Z1 + Z2 + Z3
   fit = tmle(Y = prmst, A = A, W = Xt,
              Qform = Qmodt,
              gform = Amod)
   tml = fit$estimates$ATE$psi
   
   
   
   test1 = rbind(test1, c(dr, endr, tml))
}
