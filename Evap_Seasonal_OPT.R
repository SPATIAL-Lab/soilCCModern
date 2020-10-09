## Optimizing the model for seasonal precipitation and evaporation

sm_optimizer_evap = function(Pa, Tma, PfPCQ, TmPCQ_min_a, Ra, TmWQ_min_a, DQ, spre, esw){

  #Porosity and tortuosity means of 0.35 and 0.7, respectively, and sd of 0.1.
  pores_m = 0.35
  pores_var = 0.1 ^ 2
  size = pores_m * (1 - pores_m) / pores_var - 1
  alpha = pores_m * size
  beta = (1 - pores_m) * size
  pores = rbeta(nsynth, alpha, beta)
  
  tort_m = 0.7
  tort_var = 0.1 ^ 2
  size = tort_m * (1 - tort_m) / tort_var - 1
  alpha = tort_m * size
  beta = (1 - tort_m) * size
  tort = rbeta(nsynth, alpha, beta)
  
  # Solar radiation estimate - Hargreaves and Samani (1982)
  Rs <- Ra * 0.16 * sqrt(12)
  
  #Isotope ratio constants
  RO.vsmow = 0.0020052
  RO.vpdb = 0.002067181
  
  #Climate Parameters
  PPCQ <- Pa * PfPCQ
  #At least 1mm rain in the carbonate precipitation quarter
  PPCQ = max(PPCQ, 1)
  TmPCQ <- Tma + TmPCQ_min_a
  TmPCQ_K <- TmPCQ + 273.15
  TmOOS <- (Tma * 4 -  TmPCQ) / 3
  
  #Depth to carbonate
  #top of Bk equation based on Retallack 2005 data
  z_min <- Pa * 0.0925 + 13.4
  #thickness of Bk, CMP as proxy for seasonality 
  z_thick = abs(PPCQ - Pa / 4) * 0.74 + 17.3 
  #find middle of Bk in cm
  z_mean = z_min + z_thick / 2 
  #gamma scale parameter, using 22cm as variance from Retallack fit
  theta = 22 ^ 2 / z_mean  
  #gamma shape parameter
  shp = z_mean / theta
  z = rgamma(nsynth, shape = shp, scale = theta)
  #depth in meters
  z_m <- z / 100
  
  #Soil temperatures at depth z
  #Assume thermal conductivity = 0.0007 cal / cm2  s  *C, volumetric heat capacity of 0.3 cal / cm2 *C, Quade 2013
  d <- sqrt((2 * 0.0007) / ((2 * 3.1415 / 3.154e7) * 0.3))
  #T is highest (avg summer temps) at t = 0.29, mean at t = 0 (avg spring and fall temps), low at t = 0.79 (avg winter temps)
  t <- ifelse(DQ == 3, 0.29, ifelse(DQ == 1, 0.79, 0))
  TmPCQ_soil <- Tma + (TmWQ_min_a * sin((2*3.1415) * t - z / d) / exp(z / d)) 
  TmPCQ_soil_K <- TmPCQ_soil + 273.15
  
  #Relative humidity based on quarterly precip
  h_m <- 0.25 + 0.7 * (PPCQ / 900)
  #Large variance
  h_var = 0.1 ^ 2
  size = h_m * (1 - h_m) / h_var - 1
  alpha = h_m * size
  beta = (1 - h_m) * size
  h = rbeta(nsynth, alpha, beta)
  RH <- h * 100
  
  #Potential Evapotranspiration
  #PET in mm/day, Turc 1961
  #Annual
  PET_A_D_m <- ifelse (RH < 50, 0.013 * (Tma / (Tma + 15)) * (23.8856 * Rs + 50)* (1 + ((50 - RH) / 70)), 0.0133 * (Tma / (Tma + 15)) * (23.885 * Rs + 50))
  PET_A_D_m <- pmax(PET_A_D_m, 0.01)
  theta = 0.26 ^ 2 / PET_A_D_m  
  shp = PET_A_D_m / theta
  PET_A_D = rgamma(nsynth, shape = shp, scale = theta)
  PET_A_A <- PET_A_D * 365
  
  #PCQ
  PET_PCQ_D_m <- ifelse (RH < 50, 0.013 * (TmPCQ / (TmPCQ + 15)) * (23.8856 * Rs + 50)* (1 + ((50 - RH) / 70)), 0.0133 * (TmPCQ / (TmPCQ + 15)) * (23.885 * Rs + 50))
  PET_PCQ_D_m <- pmax(PET_PCQ_D_m, 0.01)
  theta = 0.26 ^ 2 / PET_PCQ_D_m  
  shp = PET_PCQ_D_m / theta
  PET_PCQ_D = rgamma(nsynth, shape = shp, scale = theta)
  PET_PCQ <- PET_PCQ_D * 90
  
  # Calculating L (production depth) calculated from aridity index Yang et al. (2016) - approxiTmaed from Figure 4a. Assumes production depth is equal to rooting depth - Quade et al. (2007).
  AI = PET_A_A / Pa
  L = ifelse(AI > 1.4, 60, -200 * AI ^ 2 + 250 * AI + 100)
  
  # Production depth relationship with mean production depth and characteristic production depth - Quade et al. (2007)
  k = L / 2 / log(2)
  
  #Actual Evapotranspiration
  #This noise parmeter limits ETA<CMP_mm but allows variation around PET, as observed
  AET_var = rnorm(nsynth, 1, 0.2) 
  #AET in mm/quarter from Budyko curve
  AET_PCQ = PPCQ * (1 / (sqrt(1 + (1 / ((PET_PCQ / (PPCQ)) * AET_var)) ^ 2))) 
  
  #Free air porosity
  #Scales volumetrically w/ excess precipitation relative to pore space, assume a minimum of 5% volumetric water content
  FAP <- pmin((pores - ((PPCQ - AET_PCQ) / (L * 10 * pores))), pores - 0.05)
  #At least 1% free air porosity
  FAP = pmax(FAP, 0.01)
  
  #Precipitation O isotope ratios - WIDB, waterisotopes.org 
  PfPCQ = max(0.0001, PfPCQ)
  dO_P_OOS <- -15.57 + 0.60 * TmOOS
  dO_P_PCQ <- -15.57 + 0.60 * TmPCQ
  dO_P_soil_m <- (dO_P_PCQ * PfPCQ + dO_P_OOS * (1 - spre) * (1 - PfPCQ)) / (PfPCQ + (1 - spre) * (1 - PfPCQ))
  dO_P_soil = rnorm(nsynth, dO_P_soil_m, 3.04)
  R_O_P_soil = (dO_P_soil / 1000 + 1) * RO.vsmow
  R_O_P_PCQ = (dO_P_PCQ / 1000 + 1) * RO.vsmow
  
  #Atmospheric vapor O isotope ratios
  A_atmH2O_P <- 2.71828 ^ ((5.9702e6 / TmPCQ_K ^ 2 - 3.2801e4 / TmPCQ_K + 52.227) / 1000)
  R_O_atmH2O <- R_O_P_PCQ / A_atmH2O_P
  dO_atmH2O_m <- (R_O_atmH2O / RO.vsmow - 1) * 1000
  dO_atmH2O = rnorm(nsynth, dO_atmH2O_m, 1)
  
  #Soil water evaporation
  #evap is 6%+/-4% of total ET globally Good et. al. 2010
  e_mean = 0.06  
  e_var = 0.04^2
  size = e_mean*(1-e_mean)/e_var - 1
  alpha = e_mean * size
  beta = (1-e_mean) * size
  E = rbeta(nsynth, alpha, beta) * AET_PCQ
  #minimum of 1 mm of evap / quarter
  E = pmax(E, 1)
  
  #Soil water diffusion evaporation balance
  #evaporation in m/sec
  E_s <- E / (1000 * 30 * 24 * 3600) 
  #Diffusion, scaled to temperature, soil water content and tortuosity
  DIFO <- 1.637e-8 * (TmPCQ_soil_K / 216.25 - 1) ^ 2.074 * (pores - FAP) * tort
  #mean penetration depth of evap, in m
  z_i <- DIFO / E_s 
  #Diffusion ratio factor
  DRF <- 1 + 0.8 * (1 / 0.9723 - 1)
  #Surface water isotopes
  R_O_surface <- ((1 - h) * DRF * R_O_P_soil + h * R_O_atmH2O) / (1 / A_atmH2O_P)
  #Water isotopes at depth
  R_O_soil <- ((R_O_surface - R_O_P_soil) * 2.71828 ^ (-z_m / z_i)) + R_O_P_soil
  #soil water is esw % evaporated fraction
  R_O_soil = R_O_soil * esw + R_O_P_soil * (1 - esw)
  dO_soil <- (R_O_soil/RO.vsmow - 1) * 1000
    
  #Soil carbonate O isotope fractionation - Kim and O'neal 1997
  A_H20_Carb <- 2.71828 ^ ((1.803e4 / TmPCQ_soil_K - 32.42) / 1000)
  R_O_Carb <- R_O_soil * A_H20_Carb
  
  #Carb isotope values
  dO_Carb <- (R_O_Carb / RO.vpdb - 1) * 1000
  
  return(median(dO_Carb))
}

#Set up evap optimization
s = seq(0,1,0.05)
ss = seq(0, 1.05 - 0.05/21, 0.05 / 21)
ss = ss*20
ss = trunc(ss)
ss = ss/20

#Warm quarter evap optimization

parms_wq_O = data.frame(spres = ss, esws = rep(s, 21), rmse = numeric(441))

for(j in 1:nrow(parms_wq_O)){
  opt = numeric()
  for(i in 1: nrow(sites)){
    opt[i] = sm_optimizer_evap(sites$Pa[i], sites$Tma[i], sites$PfWQ[i], sites$TmWQ_min_a[i], sites$Ra[i], sites$TmWQ_min_a[i], 3, parms_wq_O$spres[j], parms_wq_O$esws[j])
  }
  
  opt = data.frame(Site = sites$Site, d18O = opt)
  opt = merge.data.frame(opt, data.aves, by.x = "Site", by.y = "Site", all.x=TRUE)
  
  mse = (opt$d18O - opt$d18O.measured)^2
  mse = mean(mse)
  parms_wq_O$rmse[j] = sqrt(mse)
}

min(parms_wq_O$rmse)
rmses_wq = matrix(parms_wq_O$rmse, 21, 21)
rmses_wq = rmses_wq[c(21:1),]
rmses.rast_wq = raster(rmses_wq, xmn=0, xmx=1, ymn=0, ymx=1)
spre_opt_wq = parms_wq_O[which.min(parms_wq_O$rmse),"spres"]
esw_opt_wq = parms_wq_O[which.min(parms_wq_O$rmse),"esws"]

#Dry quarter evap optimization

parms_dq_O = data.frame(spres = ss, esws = rep(s, 21), rmse = numeric(441))

for(j in 1:nrow(parms_dq_O)){
  opt = numeric()
  for(i in 1: nrow(sites)){
    opt[i] = sm_optimizer_evap(sites$Pa[i], sites$Tma[i], sites$PfDQ[i], sites$TmDQ_min_a[i], sites$Ra[i], sites$TmWQ_min_a[i], sites$DQ[i], parms_dq_O$spres[j], parms_dq_O$esws[j])
  }
  
  opt = data.frame(Site = sites$Site, d18O = opt)
  opt = merge.data.frame(opt, data.aves, by.x = "Site", by.y = "Site", all.x=TRUE)
  
  mse = (opt$d18O - opt$d18O.measured)^2
  mse = mean(mse)
  parms_dq_O$rmse[j] = sqrt(mse)
}

spre_opt_dq = parms_dq_O[which.min(parms_dq_O$rmse),"spres"]
esw_opt_dq = parms_dq_O[which.min(parms_dq_O$rmse),"esws"]
rmses_dq = matrix(parms_dq_O$rmse, 21, 21)
rmsess_dq = rmses_dq[c(21:1),]
rmses.rast_dq = raster(rmses_dq, xmn=0, xmx=1, ymn=0, ymx=1)
View(parms_wq_O)
spre_opt_wq
esw_opt_wq
spre_opt_dq
esw_opt_dq

#Graph

#jpeg("Evap_OPT.jpg", units="in", res=300, width=9.15, height=4.68)

layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(1.001,1))

par(mar=c(5.5, 5, 2, 0))
plot(rmses.rast_wq, main="", cex.lab=1, cex.axis=1.3, xlab="SP (Seasonal Precipitation Bias, Fraction)", ylab=expression(paste("f"[evap]," (Evaporated Water Fraction)")), xlim=c(0,1), ylim=c(0,1), zlim=c(3.5,6), legend = F)
mtext("a",3,line=0.25,adj=0, cex=1.3)
par(mar=c(5.5, 0, 2, 7))
plot(rmses.rast_dq, yaxt='n', cex.lab=1, cex.axis=1.3, main="", xlab="SP (Seasonal Precipitation Bias, Fraction)", ylab="", xlim=c(0,1), ylim=c(0,1), zlim=c(3.5,6))  #now need to make a nice plot...
mtext("b",3,line=0.25,adj=0, cex=1.3)
mtext("RMSE, (\u2030)", 4, line=0.5, cex=1.3)

#dev.off()



