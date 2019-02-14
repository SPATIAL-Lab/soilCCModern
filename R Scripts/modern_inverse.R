library(R2OpenBUGS)
library(coda)
library(rjags)
library(R2jags)

source("soil_model.R")

####analysis

parameters <- c("pCO2", "deltaA", "MAT", "T_seas", "MAP", "P_seas", "R_month", "ST", "z")
set.seed(51234)

#Lets run the bunch!
#This uses data read into the data.comp DF via 180918_forward_validation.R

posts = list()

for(i in 1: nrow(data.comp)){
  carbonate_data = list(dC_Carb=data.comp$d13C..measured.[i], dO_Carb=data.comp$d18O..measured.[i], stdevC = 0.5, stdevO = 0.5)

  bayes_fit <- jags(model.file = textConnection(soil_model), parameters.to.save = parameters, 
                    data = carbonate_data, inits = NULL, 
                    n.chains=3, n.iter = 100000, n.burnin = 10000, n.thin = 25)
  
  post = as.data.frame(bayes_fit$BUGSoutput$sims.list)
  
  posts[[i]] = post
}

ptil = numeric()
ttil = numeric()
psil = numeric()
tsil = numeric()
len = length(posts[[i]]$MAP)
for(i in 1:nrow(data.comp)){
  ptil[i] = length(posts[[i]]$MAP[posts[[i]]$MAP<data.comp$map.wc[i]])/len
  ttil[i] = length(posts[[i]]$MAT[posts[[i]]$MAT<data.comp$mat.wc[i]])/len
  psil[i] = length(posts[[i]]$MAP[posts[[i]]$P_seas<data.comp$hqp.frac[i]])/len
  tsil[i] = length(posts[[i]]$MAP[posts[[i]]$T_seas<data.comp$hqt.offset[i]])/len
}

#####Summaries and stats

med.t = numeric()
low.t = numeric()
high.t = numeric()
for(i in 1:nrow(data.comp)){
  med.t[i] = median(posts[[i]]$MAT)
  low.t[i] = quantile(posts[[i]]$MAT, probs=c(0.1))
  high.t[i] = quantile(posts[[i]]$MAT, probs=c(0.9))
}
data.comp$med.t=med.t
data.comp$low.t=low.t
data.comp$high.t=high.t
summary(lm(med.t~data.comp$mat.wc))
la.t = high.t-data.comp$mat.wc
di.t = data.comp$mat.wc-low.t
ladi.t = la.t*di.t
length(ladi.t[ladi.t>0])/length(ladi.t)

med.p = numeric()
low.p = numeric()
high.p = numeric()
for(i in 1:nrow(data.comp)){
  med.p[i] = median(posts[[i]]$MAP)
  low.p[i] = quantile(posts[[i]]$MAP, probs=c(0.1))
  high.p[i] = quantile(posts[[i]]$MAP, probs=c(0.9))
}
data.comp$med.p=med.p
data.comp$low.p=low.p
data.comp$high.p=high.p
summary(lm(med.p~data.comp$map.wc))
la.p = high.p-data.comp$map.wc
di.p = data.comp$map.wc-low.p
ladi.p = la.p*di.p
length(ladi.p[ladi.p>0])/length(ladi.p)

med.ts = numeric()
low.ts = numeric()
high.ts = numeric()
for(i in 1:nrow(data.comp)){
  med.ts[i] = median(posts[[i]]$T_seas+posts[[i]]$MAT)
  low.ts[i] = quantile(posts[[i]]$T_seas+posts[[i]]$MAT, probs=c(0.1))
  high.ts[i] = quantile(posts[[i]]$T_seas+posts[[i]]$MAT, probs=c(0.9))
}
data.comp$med.ts=med.ts
data.comp$low.ts=low.ts
data.comp$high.ts=high.ts
ts = data.comp$mat.wc+data.comp$hqt.offset
summary(lm(med.t~ts))
la.ts = high.ts-ts
di.ts = ts-low.ts
ladi.ts = la.ts*di.ts
length(ladi.ts[ladi.ts>0])/length(ladi.ts)

med.ps = numeric()
low.ps = numeric()
high.ps = numeric()
for(i in 1:nrow(data.comp)){
  med.ps[i] = median(posts[[i]]$MAP * posts[[i]]$P_seas)
  low.ps[i] = quantile(posts[[i]]$MAP * posts[[i]]$P_seas, probs=c(0.1))
  high.ps[i] = quantile(posts[[i]]$MAP * posts[[i]]$P_seas, probs=c(0.9))
}
data.comp$med.ps=med.ps
data.comp$low.ps=low.ps
data.comp$high.ps=high.ps
ps = data.comp$map.wc * data.comp$dqp.frac
summary(lm(med.ps~ps))
la.ps = high.ps-ps
di.ps = ps-low.ps
ladi.ps = la.ps*di.ps
length(ladi.ps[ladi.ps>0])/length(ladi.ps)

#####Plotting

jpeg("inverse_val.jpg", res=300, units="in", width=10, height=5)
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
par(mar=c(5,5,1,1))

plot(data.comp$mat.wc, med.t, ylim=c(-8,28), pch=16, xlab="Observed MAT", ylab="Estimated MAT")
arrows(data.comp$mat.wc, med.t, data.comp$mat.wc, high.t, angle=90, length=0.03)
arrows(data.comp$mat.wc, med.t, data.comp$mat.wc, low.t, angle=90, length=0.03)
abline(0,1)
text(-4, 27, "A")

plot(data.comp$map.wc, med.p, ylim=c(0,900), pch=16, xlab="Observed MAP", ylab="Estimated MAP")
arrows(data.comp$map.wc, med.p, data.comp$map.wc, high.p, angle=90, length=0.03)
arrows(data.comp$map.wc, med.p, data.comp$map.wc, low.p, angle=90, length=0.03)
abline(0,1)
text(80, 870, "B")

dev.off()

plot(ts, med.ts, ylim=c(0,35), pch=16, xlab="Observed CQT", ylab="Estimated CQT")
arrows(ts, med.ts, ts, high.ts, angle=90, length=0.03)
arrows(ts, med.ts, ts, low.ts, angle=90, length=0.03)
abline(0,1)

plot(ps, med.ps, ylim=c(0,200), pch=16, xlab="Observed CQP", ylab="Estimated CQP")
arrows(ps, med.ps, ps, high.ps, angle=90, length=0.03)
arrows(ps, med.ps, ps, low.ps, angle=90, length=0.03)
abline(0,1)


#Now we'll do it with clumped!
#This uses data read into the data.comp DF via 180918_forward_validation.R

source("soil_model_clumped.R")
parameters <- c("pCO2", "deltaA", "MAT", "T_seas", "MAP", "P_seas", "R_month", "ST", "z")
set.seed(21334)

posts.cl = list()

for(i in 1: nrow(data.clump)){
  carbonate_data = list(dC_Carb=data.clump$d13C.measured[i], 
                        dO_Carb=data.clump$d18O.measured[i], 
                        D47_Carb=data.clump$D47.measured[i],
                        stdevC = 0.5, stdevO = 0.5, stdev47 = 0.0012)
  
  bayes_fit <- jags(model.file = textConnection(soil_model_clumped), parameters.to.save = parameters, 
                    data = carbonate_data, inits = NULL, 
                    n.chains=3, n.iter = 100000, n.burnin = 10000, n.thin = 25)
  
  post = as.data.frame(bayes_fit$BUGSoutput$sims.list)
  posts.cl[[i]] = post
}

#####Summaries and stats

med.t = numeric()
low.t = numeric()
high.t = numeric()
for(i in 1:length(posts.cl)){
  med.t[i] = median(posts.cl[[i]]$MAT)
  low.t[i] = quantile(posts.cl[[i]]$MAT, probs=c(0.1))
  high.t[i] = quantile(posts.cl[[i]]$MAT, probs=c(0.9))
}
data.clump$med.t=med.t
data.clump$low.t=low.t
data.clump$high.t=high.t
summary(lm(med.t~data.clump$mat.wc))
la.t = high.t-data.clump$mat.wc
di.t = data.clump$mat.wc-low.t
ladi.t = la.t*di.t
length(ladi.t[ladi.t>0])/length(ladi.t)

med.p = numeric()
low.p = numeric()
high.p = numeric()
for(i in 1:length(posts.cl)){
  med.p[i] = median(posts.cl[[i]]$MAP)
  low.p[i] = quantile(posts.cl[[i]]$MAP, probs=c(0.1))
  high.p[i] = quantile(posts.cl[[i]]$MAP, probs=c(0.9))
}
data.clump$med.p=med.p
data.clump$low.p=low.p
data.clump$high.p=high.p
summary(lm(med.p~data.clump$map.wc))
la.p = high.p-data.clump$map.wc
di.p = data.clump$map.wc-low.p
ladi.p = la.p*di.p
length(ladi.p[ladi.p>0])/length(ladi.p)

med.ts = numeric()
low.ts = numeric()
high.ts = numeric()
for(i in 1:length(posts.cl)){
  med.ts[i] = median(posts.cl[[i]]$T_seas+posts.cl[[i]]$MAT)
  low.ts[i] = quantile(posts.cl[[i]]$T_seas+posts.cl[[i]]$MAT, probs=c(0.1))
  high.ts[i] = quantile(posts.cl[[i]]$T_seas+posts.cl[[i]]$MAT, probs=c(0.9))
}
data.clump$med.ts=med.ts
data.clump$low.ts=low.ts
data.clump$high.ts=high.ts
ts = data.clump$mat.wc+data.clump$hqt.offset
summary(lm(med.t~ts))
la.ts = high.ts-ts
di.ts = ts-low.ts
ladi.ts = la.ts*di.ts
length(ladi.ts[ladi.ts>0])/length(ladi.ts)

med.ps = numeric()
low.ps = numeric()
high.ps = numeric()
for(i in 1:length(posts.cl)){
  med.ps[i] = median(posts.cl[[i]]$MAP * posts.cl[[i]]$P_seas)
  low.ps[i] = quantile(posts.cl[[i]]$MAP * posts.cl[[i]]$P_seas, probs=c(0.1))
  high.ps[i] = quantile(posts.cl[[i]]$MAP * posts.cl[[i]]$P_seas, probs=c(0.9))
}
data.clump$med.ps=med.ps
data.clump$low.ps=low.ps
data.clump$high.ps=high.ps
ps = data.clump$map.wc * data.clump$dqp.frac
summary(lm(med.ps~ps))
la.ps = high.ps-ps
di.ps = ps-low.ps
ladi.ps = la.ps*di.ps
length(ladi.ps[ladi.ps>0])/length(ladi.ps)

######Plotting!

jpeg("inverse_val_clump.jpg", res=300, units="in", width=10, height=5)
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
par(mar=c(5,5,1,1))

plot(data.clump$mat.wc, med.t, ylim=c(-8,28), pch=16, xlab="Observed MAT", ylab="Estimated MAT")
arrows(data.clump$mat.wc, med.t, data.clump$mat.wc, high.t, angle=90, length=0.03)
arrows(data.clump$mat.wc, med.t, data.clump$mat.wc, low.t, angle=90, length=0.03)
abline(0,1)
text(-4, 27, "A")

plot(data.clump$map.wc, med.p, ylim=c(0,900), pch=16, xlab="Observed MAP", ylab="Estimated MAP")
arrows(data.clump$map.wc, med.p, data.clump$map.wc, high.p, angle=90, length=0.03)
arrows(data.clump$map.wc, med.p, data.clump$map.wc, low.p, angle=90, length=0.03)
abline(0,1)
text(80, 870, "B")

dev.off()

plot(ts, med.ts, ylim=c(0,35), pch=16, xlab="Observed CQT", ylab="Estimated CQT")
arrows(ts, med.ts, ts, high.ts, angle=90, length=0.03)
arrows(ts, med.ts, ts, low.ts, angle=90, length=0.03)
abline(0,1)

plot(ps, med.ps, ylim=c(0,200), pch=16, xlab="Observed CQP", ylab="Estimated CQP")
arrows(ps, med.ps, ps, high.ps, angle=90, length=0.03)
arrows(ps, med.ps, ps, low.ps, angle=90, length=0.03)
abline(0,1)

ce = merge.data.frame(data.comp, data.clump, by.x="Site", by.y="Site")

attach(ce)
plot(med.t.x, med.t.y, pch=16, xlim=c(-10,30), ylim=c(0,40), xlab="MAT, no clumped", ylab="MAT, with clumped")
arrows(med.t.x, med.t.y, med.t.x, low.t.y, angle=90, length=0.03)
arrows(med.t.x, med.t.y, med.t.x, high.t.y, angle=90, length=0.03)
arrows(med.t.x, med.t.y, low.t.x, med.t.y, angle=90, length=0.03)
arrows(med.t.x, med.t.y, high.t.x, med.t.y, angle=90, length=0.03)
abline(0,1)
plot(med.p.x, med.p.y, pch=16, xlim=c(0,900), ylim=c(0,800), xlab="MAP, no clumped", ylab="MAP, with clumped")
arrows(med.p.x, med.p.y, med.p.x, low.p.y, angle=90, length=0.03)
arrows(med.p.x, med.p.y, med.p.x, high.p.y, angle=90, length=0.03)
arrows(med.p.x, med.p.y, low.p.x, med.p.y, angle=90, length=0.03)
arrows(med.p.x, med.p.y, high.p.x, med.p.y, angle=90, length=0.03)
abline(0,1)
plot(med.ts.x, med.ts.y, pch=16, xlim=c(0,40), ylim=c(10,60), xlab="CQT, no clumped", ylab="CQT, with clumped")
arrows(med.ts.x, med.ts.y, med.ts.x, low.ts.y, angle=90, length=0.03)
arrows(med.ts.x, med.ts.y, med.ts.x, high.ts.y, angle=90, length=0.03)
arrows(med.ts.x, med.ts.y, low.ts.x, med.ts.y, angle=90, length=0.03)
arrows(med.ts.x, med.ts.y, high.ts.x, med.ts.y, angle=90, length=0.03)
abline(0,1)
plot(med.ps.x, med.ps.y, pch=16, xlim=c(0,200), ylim=c(0,150), xlab="CQP, no clumped", ylab="CQP, with clumped")
arrows(med.ps.x, med.ps.y, med.ps.x, low.ps.y, angle=90, length=0.03)
arrows(med.ps.x, med.ps.y, med.ps.x, high.ps.y, angle=90, length=0.03)
arrows(med.ps.x, med.ps.y, low.ps.x, med.ps.y, angle=90, length=0.03)
arrows(med.ps.x, med.ps.y, high.ps.x, med.ps.y, angle=90, length=0.03)
abline(0,1)
detach(ce)
