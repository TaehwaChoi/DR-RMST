pkg = c("survival", "pseudo", "tmle", "grf", "randomForest","tictoc", "gbm")
sapply(pkg, library, character.only=T)
#####
# test21: L=10
# test22 : L=20
simdata2 = function(n, tau) 
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
   A = rbinom(n, 1, prob = expit(x1 - x4*(x3+x5))) 
   # True parameters 
   model = (A==1)*exp(-5 - sin(2.1*x1) + 
                         0.5*x3*x4 + 1.4*x4*x5) +
      (A==0)*exp(-2 - sin(4*x1) - 1.6*x2*x4 
                 + 1.2*x3*x5)
   Time  = rexp(n, rate = model)              # (unknown) true time 
   # Given information (our dataset)   
   Cens  = rexp(n, tau)                       # Censoring rate ~ 30%         
   Y     = pmin(Time, Cens)                   # Observed time (right censored)
   delta = ifelse(Time <= Cens, 1, 0)         # mean(I(T<=C)) = Prob(T<=C)            
   range(Y[A==0]) ; range(Y[A==1])
   1-mean(delta)
   dat = data.frame(Y = round(Y,5), delta, A, x1, x2, x3, x4, x5)
}

n     = 600         # sample size
nsim  = 1000      # number of iteration
tau  = exp(-4.9)    # determines censoring rate
truepar = c(4.194, 10.728)


test21 = NULL
L     = 10        
Qmod = "prmst ~ x1+x2+x3+x4+x5"
Qmodt = "Y ~ x1+x2+x3+x4+x5"
Amod = "as.factor(A) ~ x1+x2+x3+x4+x5"
for (i in 1:nsim)
{
   set.seed(i)
   dt = simdata2(n, tau)
   Y = dt$Y; delta = dt$delta; A = dt$A; 
   Z1 = dt$Z1; Z2 = dt$Z2; Z3 = dt$Z3
   Xt = with(dt, as.matrix(data.frame(x1,x2,x3,x4,x5)))
   Xc = with(dt, as.matrix(data.frame(x1,x2,x3,x4,x5)))
   
   dt$prmst = prmst = pseudomean(time = Y, event = delta, tmax = L)
   m1 = predict(lm(Qmod, data = subset(dt, A == 1)), newdata = dt)
   m0 = predict(lm(Qmod, data = subset(dt, A == 0)), newdata = dt)
   # Reg-L
   reg = mean(m1-m0)
   
   # dr-L
   pi1 = fitted(glm(Amod, data = dt, family = binomial))
   dr = mean((A/pi1 - (1-A)/(1-pi1))*prmst - 
                (A-pi1)/pi1*m1 - (A-pi1)/(1-pi1)*m0)
   
   
   # dr+-L
   idx1 = which(A==1); idx0 = which(A==0)
   d = log(mean(A*(1+exp(pi1)) / exp(pi1)))
   ENpi1 = pi1
   wt1 = ((1-ENpi1)/ENpi1^2)[idx1]
   wt0 = ((ENpi1)/(1-ENpi1)^2)[idx0]
   EN1 = predict(lm(Qmod, data = subset(dt, A ==1), weights = wt1), dt)
   EN0 = predict(lm(Qmod, data = subset(dt, A ==0), weights = wt0), dt)
   endr = mean((A/ENpi1-(1-A)/(1-ENpi1))*prmst -
                  (A-ENpi1)/ENpi1*EN1 - (A-ENpi1)/(1-ENpi1)*EN0)
   
   # ipw-gbm
   gfit = gbm((A) ~ x1+x2+x3+x4+x5, data = dt, distribution = "bernoulli")
   pigb = predict(gfit, dt, n.trees = 100, type = "response")
   ipgb = mean((A/pigb - (1-A)/(1-pigb))*prmst)
   
   
   # dr-gbm
   drgb = mean((A/pigb - (1-A)/(1-pigb))*prmst - 
                  (A-pigb)/pigb*m1 - (A-pigb)/(1-pigb)*m0)
   
   
   # TMLE-gbm
   fit = tmle(Y = prmst, A = A, W = Xt,
              Qform = Qmodt,
              gform = Amod,
              g.SL.library = "SL.gbm")
   tmgb = fit$estimates$ATE$psi
   
   
   #TMLE-sl
   fit = tmle(Y = prmst, A = A, W = Xt,
              Q.SL.library = c("SL.glm", "SL.glm.interaction", "SL.step"),
              g.SL.library =  c("SL.glm", "SL.glm.interaction", "SL.step"))
   tmsl = fit$estimates$ATE$psi
   
   # CF
   forest = causal_forest(X = Xc, Y = prmst, W = A)
   cf =mean(predict(forest)$pred)
   
   
   test21 = rbind(test21, c(reg, dr, endr, ipgb, drgb, tmgb, tmsl, cf))
}



test22 = NULL
L     = 20        
Qmod = "prmst ~ x1+x2+x3+x4+x5"
Qmodt = "Y ~ x1+x2+x3+x4+x5"
Amod = "as.factor(A) ~ x1+x2+x3+x4+x5"
for (i in 1:nsim)
{
   dt = simdata2(n, tau)
   Y = dt$Y; delta = dt$delta; A = dt$A; 
   Z1 = dt$Z1; Z2 = dt$Z2; Z3 = dt$Z3
   Xt = with(dt, as.matrix(data.frame(x1,x2,x3,x4,x5)))
   Xc = with(dt, as.matrix(data.frame(x1,x2,x3,x4,x5)))
   
   dt$prmst = prmst = pseudomean(time = Y, event = delta, tmax = L)
   m1 = predict(lm(Qmod, data = subset(dt, A == 1)), newdata = dt)
   m0 = predict(lm(Qmod, data = subset(dt, A == 0)), newdata = dt)
   # Reg-L
   reg = mean(m1-m0)
   
   # dr-L
   pi1 = fitted(glm(Amod, data = dt, family = binomial))
   dr = mean((A/pi1 - (1-A)/(1-pi1))*prmst - 
                (A-pi1)/pi1*m1 - (A-pi1)/(1-pi1)*m0)
   
   
   # dr+-L
   idx1 = which(A==1); idx0 = which(A==0)
   d = log(mean(A*(1+exp(pi1)) / exp(pi1)))
   ENpi1 = pi1
   wt1 = ((1-ENpi1)/ENpi1^2)[idx1]
   wt0 = ((ENpi1)/(1-ENpi1)^2)[idx0]
   EN1 = predict(lm(Qmod, data = subset(dt, A ==1), weights = wt1), dt)
   EN0 = predict(lm(Qmod, data = subset(dt, A ==0), weights = wt0), dt)
   endr = mean((A/ENpi1-(1-A)/(1-ENpi1))*prmst -
                  (A-ENpi1)/ENpi1*EN1 - (A-ENpi1)/(1-ENpi1)*EN0)
   
   # ipw-gbm
   gfit = gbm((A) ~ x1+x2+x3+x4+x5, data = dt, distribution = "bernoulli")
   pigb = predict(gfit, dt, n.trees = 100, type = "response")
   ipgb = mean((A/pigb - (1-A)/(1-pigb))*prmst)
   
   
   # dr-gbm
   drgb = mean((A/pigb - (1-A)/(1-pigb))*prmst - 
                  (A-pigb)/pigb*m1 - (A-pigb)/(1-pigb)*m0)
   
   
   # TMLE-gbm
   fit = tmle(Y = prmst, A = A, W = Xt,
              Qform = Qmodt,
              gform = Amod,
              g.SL.library = "SL.gbm")
   tmgb = fit$estimates$ATE$psi
   
   
   #TMLE-sl
   fit = tmle(Y = prmst, A = A, W = Xt,
              Q.SL.library = c("SL.glm", "SL.glm.interaction", "SL.step"),
              g.SL.library =  c("SL.glm", "SL.glm.interaction", "SL.step"))
   tmsl = fit$estimates$ATE$psi
   
   # CF
   forest = causal_forest(X = Xc, Y = prmst, W = A)
   cf =mean(predict(forest)$pred)
   
   
   test22 = rbind(test22, c(reg, dr, endr, ipgb, drgb, tmgb, tmsl, cf))
}

