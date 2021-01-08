setwd(getwd())

rm(list=ls())
library(pseudo)
library(dplyr)

load("nhefs.RData")
range(nhefs$yrdth) # 1983 - 1992
nhefs$survtime = (nhefs$yrdth - 83)*12 + nhefs$modth # survtime month로 단위 바꾸기  
range(nhefs$survtime) 
sum(is.na(nhefs$survtime)) #0

nhefs = within(nhefs, {
   age = as.numeric(age)
   education = as.factor(education)
   active = as.factor(active)
   exercise = as.factor(exercise)
})

### 
library(tmle)
library(SuperLearner)
library(grf)
library(gbm)
L1 = 12*5  # administrative censoring at 60 months 
L2 = 12*10 # administrative censoring at 120 months
quantile(sort(unique(nhefs$survtime)), probs = seq(0,1,0.2))
Lvec1  = c(1, seq(12, 120, by = 12))     # Lvec1 
Lvec2  = c(1, seq(12, 60, by = 12))     # Lvec2 

# (A, Y, C, X, Pseudo)
A = nhefs$qsmk
Y = nhefs$survtime
death = nhefs$death
n = nrow(nhefs)
dt = nhefs
dt$prmst = pseudomean(time = Y, event = death, tmax = L2)

set.seed(123)
gfit = gbm(qsmk ~ sex + race + age + as.factor(education)
           + smokeintensity + smokeyrs + as.factor(exercise)
           + as.factor(active) + wt71, data = dt, distribution = "bernoulli")
pigb = predict(gfit, dt, n.trees = 100, type = "response")

mod.m0 = predict(gbm(prmst ~ sex + race + age + as.factor(education)
                     + smokeintensity + smokeyrs + as.factor(exercise)
                     + as.factor(active) + wt71, data = subset(dt, A==0), distribution = "gaussian"), dt, n.trees = 100)
mod.m1 = predict(gbm(prmst ~ sex + race + age + as.factor(education)
                     + smokeintensity + smokeyrs + as.factor(exercise)
                     + as.factor(active) + wt71, data = subset(dt, A==1), distribution = "gaussian"), dt, n.trees = 100)

(dr3 = mean(with(dt, (A/pigb - (1-A)/(1-pigb))*prmst - (A-pigb)/pigb*mod.m1 -
                    (A-pigb)/(1-pigb)*mod.m0)))


Xt = model.matrix(~-1+sex + race + age + as.factor(education)
                  + smokeintensity + smokeyrs + as.factor(exercise)
                  + as.factor(active) + wt71, data = dt)
cova = with(dt, cbind(sex, race, age, education, smokeintensity, smokeyrs, exercise,
                      active, wt71))
fit = tmle(Y = dt$prmst, A = A, W = cova,
           Q.SL.library = c("SL.glm", "SL.glm.interaction", "SL.step"),
           g.SL.library =  c("SL.glm", "SL.glm.interaction", "SL.step"))
tmsl = fit$estimates$ATE$psi


forest = causal_forest(X = Xt, Y = dt$prmst, W = A)
cf =mean(predict(forest)$pred)
adv = c(dr3, tmsl, cf)

test = NULL
for (i in 1:100) {
   set.seed(i)
   n = nrow(nhefs)
   idx = sample(1:n,n,replace=T)
   dt.boot = dt[idx,]
   dt.boot$prmst.boot = with(dt.boot, pseudomean(survtime, qsmk, L2))
   
   mod.m0 = predict(gbm(prmst.boot ~ sex + race + age + as.factor(education)
                        + smokeintensity + smokeyrs + as.factor(exercise)
                        + as.factor(active) + wt71, data = subset(dt.boot, A==0), distribution = "gaussian"),
                    dt.boot, n.trees = 100)
   mod.m1 = predict(gbm(prmst.boot ~ sex + race + age + as.factor(education)
                        + smokeintensity + smokeyrs + as.factor(exercise)
                        + as.factor(active) + wt71, data = subset(dt.boot, A==1), distribution = "gaussian"),
                    dt.boot, n.trees = 100)
   
   set.seed(i)
   gfit = gbm((qsmk) ~ sex + race + age + as.factor(education)
              + smokeintensity + smokeyrs + as.factor(exercise)
              + as.factor(active) + wt71, data = dt.boot, distribution = "bernoulli")
   pigb = predict(gfit, dt.boot, n.trees = 100, type = "response")
   dr3 = mean(with(dt.boot, (A.boot/pigb - (1-A.boot)/(1-pigb))*prmst.boot - 
                      (A.boot-pigb)/pigb*mod.m1 -
                      (A.boot-pigb)/(1-pigb)*mod.m0))
   
   
   Xt.boot = Xt[idx,]
   A.boot = A[idx]
   cova.boot = cova[idx,]
   fit = tmle(Y = dt.boot$prmst.boot, A = A.boot, W = cova.boot,
              Q.SL.library = c("SL.glm", "SL.glm.interaction", "SL.step"),
              g.SL.library =  c("SL.glm", "SL.glm.interaction", "SL.step"))
   tmsl = fit$estimates$ATE$psi
   
   
   
   forest = causal_forest(X = Xt.boot, Y = dt.boot$prmst.boot, W = A[idx])
   cf =mean(predict(forest)$pred)
   test = rbind(test, cbind(dr3, tmsl, cf))
}

xtable::xtable(cbind(
   adv,
   apply(test,2,sd),
   2*sapply(-abs(adv) /  apply(test,2,sd), function(x) pnorm(x))), digits = 3)

round(1.96*2*apply(test,2,sd), 3)

