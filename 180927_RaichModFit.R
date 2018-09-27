resp_mod = " model {


  for(i in 1:nobs){
    theta.P[i] = P.var / P.m[i]
    k.P[i] = P.m[i] / theta.P[i]    
    P[i] ~ dgamma(k.P[i], 1 / theta.P[i])
    Ta[i] ~ dnorm(Ta.m[i], 1 / T.var)
    Rs_ann.m[i] = F * exp(Q * Ta[i]) * P[i] / (K + P[i]) 
    Rs_ann[i] ~ dnorm(Rs_ann.m[i], 1 / R.var) 
  }
  
  R.var = 50 ^ 2
  T.var = 1 ^ 2
  P.var = 25 ^ 2

  F ~ dnorm(F.m, 1 / F.sd ^ 2)
  F.sd = 0.25
  F.m = 1.25

  Q ~ dnorm(Q.m, 1 / Q.sd ^ 2)
  Q.sd = 0.015
  Q.m = 0.0545

  K ~ dgamma(k.K, theta.K)
  k.K = K.m / theta.K
  theta.K = K.sd ^ 2 / K.m
  K.sd = 2
  K.m = 4.259
}
"

library(R2OpenBUGS)
library(coda)
library(rjags)
library(R2jags)

parameters = c("K", "Q", "F")
set.seed(1395)

d = read.csv("srdb-data.csv")
d = d[!is.na(d$Rs_annual),]
d = d[!is.na(d$MAP),]
d = d[!is.na(d$MAT),]
d = d[d$Ecosystem_state == "Natural"]
d = d[d$Manipulation == "None",]

rdat = list(P.m = d$MAP/120, Ta.m = d$MAT, Rs_ann = d$Rs_annual/365, nobs=nrow(d))

rmod <- jags(model.file = textConnection(resp_mod), parameters.to.save = parameters, 
                  data = rdat, inits = NULL, 
                  n.chains=3, n.iter = 5000, n.burnin = 100, n.thin = 5)

rmod
rmod.mcmc = as.mcmc(rmod)
plot(rmod.mcmc)
