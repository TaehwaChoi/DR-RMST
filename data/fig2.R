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
glimpse(nhefs)
range(nhefs$age) # age from  25 to 74 (baseline age)




## 2. Pseudo-observation approach
# 2-1 pseudo/sw.a 
# Generate psuedo values 
L1 = 12*5  # administrative censoring at 60 months 
L2 = 12*10 # administrative censoring at 120 months
quantile(sort(unique(nhefs$survtime)), probs = seq(0,1,0.2))
Lvec1  = c(1, seq(12, 120, by = 12))     # Lvec1 (11개 시점)
Lvec2  = c(1, seq(12, 60, by = 12))     # Lvec2 (11개 시점)

# (A, Y, C, X, Pseudo)
A = nhefs$qsmk
Y = nhefs$survtime
death = nhefs$death
n = nrow(nhefs)
dt = nhefs
dt$prmst = pseudomean(time = Y, event = death, tmax = L1)
#dt$prmst = pseudomean(time = Y, event = death, tmax = L2)



### survminer 
library(survminer)
library(survival)
library(ggplot2)


nhefs$qsmk = factor(nhefs$qsmk)
nfit = survfit(Surv(survtime, death) ~ as.factor(qsmk), data = nhefs)
temp = data.frame(qsmk = levels(nhefs$qsmk))
nfit = survfit(coxph(Surv(survtime, death) ~ qsmk, data = nhefs), temp)
nfit1 = as.data.frame(cbind(time = nfit$time, surv0 = nfit$surv[,1], surv1 = nfit$surv[,2]))

mod.p = glm(qsmk ~ sex + race + age + as.factor(education)
            + smokeintensity + smokeyrs + as.factor(exercise)
            + as.factor(active) + wt71, data = dt, family=binomial())
p.1 <- predict(mod.p, type = "response")  # p.1 = Prob(A=1|Z)

utime = sort(unique(nhefs$survtime))
ppp = pseudosurv(time = Y, event = death, tmax = utime)$pseudo

dr1 = dr0 = 0
for (j in 1:length(utime)) {
   dtt = cbind(ppp = ppp[,j], dt)
   m11 <- predict(lm(ppp ~ sex + race + age + as.factor(education)
                     + smokeintensity + smokeyrs + as.factor(exercise)
                     + as.factor(active) + wt71, data = subset(dtt, A==1)), dtt)
   m00 <- predict(lm(ppp ~ sex + race + age + as.factor(education)
                     + smokeintensity + smokeyrs + as.factor(exercise)
                     + as.factor(active) + wt71, data = subset(dtt, A==0)), dtt)
   dr1[j] = mean(A/p.1*ppp[,j] - (A-p.1)/p.1*m11)
   dr0[j] = mean((1-A)/(1-p.1)*ppp[,j] + ((A-p.1)/(1-p.1)*m00))
}


nfit2 = as.data.frame(cbind(time = utime, surv0 = dr0, surv1 = dr1))
a1 = ggplot(nfit1, aes(x=time,  y=surv)) + 
   geom_line(aes(y = surv0, colour = "0")) + 
   geom_line(aes(y = surv1, colour = "1")) +
   xlab("Months") + ylab("Survival Probability") + 
   scale_x_continuous(limits = c(0, 120), breaks=seq(0,120,12)) +
   scale_y_continuous(limits=c(0.6, 1), breaks=seq(0.6, 1, 0.1)) +
   scale_shape_discrete(labels = c("Smoker", "Non-smoker")) +
   scale_colour_discrete(labels = c("Smoker", "Non-smoker")) + 
   ggtitle("(a) Unadjusted KM") + 
   labs(colour="") +
   theme_bw() + 
   theme(legend.position="bottom", plot.title = element_text(size = 15)) 
a2 = ggplot(nfit2, aes(x=time,  y=surv)) + 
   geom_line(aes(y = surv0, colour = "0")) + 
   geom_line(aes(y = surv1, colour = "1")) +
   xlab("Months") + ylab("Survival Probability") + 
   scale_x_continuous(limits = c(0, 120), breaks=seq(0,120,12)) +
   scale_y_continuous(limits=c(0.6, 1), breaks=seq(0.6, 1, 0.1)) +
   scale_shape_discrete(labels = c("Smoker", "Non-smoker")) +
   scale_colour_discrete(labels = c("Smoker", "Non-smoker")) + 
   ggtitle("(b) Covariate-adjusted KM") + 
   theme(plot.title = element_text(size = 40)) +
   labs(colour="") +
   theme_bw() + 
   theme(legend.position="bottom", plot.title = element_text(size = 15)) 
library(patchwork)
gridExtra::grid.arrange(a1, a2, ncol = 2)

combined = a1 + a2 & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
