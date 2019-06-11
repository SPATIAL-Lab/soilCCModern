sm_forward_evap = function(MAP, MAT, P_seas, T_seas, Ra, pCO2){
  
  deltaA = rnorm(nsynth, -6.5, 0.3)
  #Make porosity and tortuosity scale from ST
  pores = rnorm(nsynth, 0.4, 0.1)
  tort = rnorm(nsynth, 0.7, 0.1)
  #Evaporated soil water (100%)
  esw = 1
  
  # Solar radiation estimate
  Rs <- Ra * 0.16 * sqrt(12)
  
  #Diffusion ratio of carbon isotopes 13C/12C
  DIF.ratio = 1.004443
  
  #Isotope ratio constants
  RC.vpdb = 0.011237
  RO.vsmow = 0.0020052
  RO.vpdb = 0.002067181
  
  #Basal depth
  L = 100
  
  #Climate
  CQP <- MAP * P_seas
  #At least 3mm rain in the carbonate precipitation quarter
  CQP = max(CQP, 3)
  CMP_mm <- CQP / 3
  CMP_cm <- CQP / 30
  CQT <- MAT + T_seas
  CQT_K <- CQT + 273.15
  
  #Soil carbonate C isotopes - damping temperature to soil temperature at depth z
  ## Assume thermal conductivity = 0.0007 cal / cm2  s  *C, volumetric heat capacity of 0.3 cal / cm2 *C, Quade 2013, dry sandy soils
  
  d <- sqrt((2*0.0007)/((2*3.1415/3.154e7)*0.3))
  #For sensitivity tests, t=1 - hot quarter assumed
  t <- ifelse(CQT > MAT + 3, 0.3, ifelse(CQT < MAT - 3, 0.8, 0.05))
  T_soil <- MAT + (T_seas * sin((2*3.1415) * t - z/d) / exp(z/d)) 
  T_soil_K <- T_soil + 273.15
  A_CO2_Carb <- 2.71828 ^ (-2.988e3 / T_soil_K ^ 2 + 7.6663 / T_soil_K - 0.0024612)
  
  
  #Convert pCO2 to units mol/cm^3
  pCO2_mcc = pCO2 / (0.08206 * T_soil_K * 10^9)  #mol/cm^3
  
  #Relative humidity, now beta distribution, based on quarterly precip
  h_m <- 0.25 + 0.7 * (CQP / 900)
  h_var = 0.05^2
  size = h_m*(1-h_m)/h_var - 1
  alpha = h_m * size
  beta = (1-h_m) * size
  h = rbeta(nsynth, alpha, beta)
  RH <- h * 100
  
  # Precipitation O isotope ratios - OIPC mid-latitudes
  dO_P_m <- -14.76 + 0.56 * MAT
  dO_P = rnorm(nsynth, dO_P_m, 1.46)
  R_O_P = (dO_P / 1000 + 1) * RO.vsmow
  
  #Depth to carbonate, now gamma dist
  z_min <- MAP*0.0925 + 13.4  #top of Bk equation based on Retallack 2005 data
  z_thick = abs(CMP_mm - MAP/12) * 0.74 + 17.3 #thickness of Bk, CMP as proxy for seasonality 
  z_mean = z_min + z_thick/2 #find middle of Bk in cm
  theta = 20^2/z_mean  #gamma scale parameter, using 20cm as variance from Retallack fit
  k = z_mean / theta #gamma shape parameter
  z = rgamma(nsynth, shape = k, scale = theta)
  z <- pmin(z, 100) 
  z_m <- z/100
  
  #Respiration rate, now gamma dist
  R_day_m <- 1.25 * exp(0.0545 * CQT) * CMP_cm / (4.259 + CMP_cm)  #Raich 2002, gC/m2day
  theta = (R_day_m*0.5)^2/R_day_m #gamma scale parameter, using mean residual of 50% based on Raich validation data 
  k = R_day_m / theta #gamma shape parameter
  R_day = rgamma(nsynth, shape = k, scale = theta) #lets use gamma for these quants bounded at zero....
  R_day = R_day / (12.01 * 100^2)  #molC/cm2day
  R_sec <- R_day / (24 * 3600)  #molC/cm2s
  R_sec = R_sec / L / pores #molC/cm3s
  
  #Potential ET
  ETP_D_m <- ifelse (RH < 50, 0.013 * (CQT / (CQT + 15)) * (23.8856 * Rs + 50)* (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50))
  ETP_D = rnorm(nsynth, ETP_D_m, 0.2)  #PET in mm/day, Turc 1961
  ETP_M <- ETP_D * 30  #mm/month
  
  #Actual ET
  ETA_var = rnorm(nsynth, 1, 0.2) #This noise parmeter limits ETA<CMP_mm but allows variation around ETP, as observed
  ETA = CMP_mm*3 * (1 / (sqrt(1 + (1 / ((ETP_M / (CMP_mm*3)) * ETA_var)) ^ 2))) #AET in mm/month from Budyko curve
  #here scaled eta to quarter precip, assuming potential carry-over
  
  #Free air porosity
  #Have updated, now scales volumetrically w/ excess precipitation relative to pore space, assume a minimum of 5% volumetric water content
  FAP <- pmin((pores - ((CMP_mm - ETA)/(L*10*pores))), pores-0.05)
  FAP = pmax(FAP,0.01) #dimensionless. At least 1% free air porosity
  SWC = pores - FAP
  
  #CO2 Diffusion coefficient - based on temp, FAP, tort
  DIFC = FAP * tort * 0.1369 * (T_soil_K / 273.15) ^ 1.958
  
  #Water limitation of discriminaton, Diefendorf
  W_m <- 22.65 - (1.2 * (MAP + 975)) / (27.2 + 0.04 * (MAP + 975))
  W = rnorm(nsynth, W_m, 1.1)
  
  #CO2 effect on discrimination, Schubert
  deltaP_pCO2_m <- 28.26 * 0.35 * (pCO2 + 15) / (28.26 + 0.35 * (pCO2 + 15))
  deltaP_pCO2 = rnorm(nsynth, deltaP_pCO2_m, 1.5)
  
  #Discrimination
  deltaP <- deltaA - (deltaP_pCO2 - W)
  
  #Soil CO2 C isotopes
  deltaA_hat <- (deltaA / 1000 + 1) * RC.vpdb / (1 + RC.vpdb * (deltaA / 1000 + 1))
  deltaP_hat <- (deltaP / 1000 + 1) * RC.vpdb / (1 + RC.vpdb * (deltaP / 1000 + 1))
  dC_Soil.resp = R_sec/(DIFC) * (L * z - z^2 / 2)
  dC_Soil.num = dC_Soil.resp * DIF.ratio * deltaP_hat + pCO2_mcc * deltaA_hat
  dC_Soil.denom = dC_Soil.resp * (1 - DIF.ratio * deltaP_hat) + pCO2_mcc * (1 - deltaA_hat)
  dC_Soil = (dC_Soil.num / (dC_Soil.denom * RC.vpdb) - 1) * 1000
  R_Soil <- (dC_Soil / 1000 + 1) * RC.vpdb
  R_Carb <- R_Soil / A_CO2_Carb
  
  #Soil water evaporation, now beta dist
  e_mean = 0.06  #evap is 6% of total ET globally, but could be more in arid to sub-humid climates so we estimate 10% +/- 5%
  e_var = 0.04^2
  size = e_mean*(1-e_mean)/e_var - 1
  alpha = e_mean * size
  beta = (1-e_mean) * size
  E = rbeta(nsynth, alpha, beta) * ETA #mm/month
  E = pmax(E, 1) #mm/month minimum of 1 mm of evap / month
  
  #Soil water diffusion evaporation balance
  E_s <- E / (1000 * 30 * 24 * 3600) #evaporation in m/sec
  DIFO <- 1.637e-8 * (T_soil_K / 216.25 - 1) ^ 2.074 * (pores-FAP) * tort
  z_i <- DIFO / E_s #mean penetration depth of evap, in m
  #Soil water O isotopes
  DRF <- 1 + 0.8 * (1 / 0.9723 - 1)
  R_O_surface <- ((1 - h) * DRF * R_O_P + h * R_O_atm) / (1 / A_atmP)
  R_O_soil <- ((R_O_surface - R_O_P) * 2.71828 ^ (-z_m / z_i)) + R_O_P
  R_O_soil = R_O_soil * esw + R_O_P * (1 - esw)  #soil water is esw % evaporated fraction
  dO_soil <- (R_O_soil/RO.vsmow - 1) * 1000
  
  #Soil carbonate O isotopes
  A_O <- 2.71828 ^ ((2.78e6 / T_soil_K ^ 2 - 2.89) / 1000)
  R_O_Carb <- R_O_soil * A_O
  
  dC_Carb <- (R_Carb / RC.vpdb - 1) * 1000
  dO_Carb <- (R_O_Carb / RO.vpdb - 1) * 1000
  
  dat = c(median(dC_Carb), median(dO_Carb), sd(dC_Carb), sd(dO_Carb))
  
  return(dat)
}

## w/ evap
hq_pred = data.frame(d13C=numeric(0), d18O=numeric(0), d13C_sd=numeric(0), d18O_sd=numeric(0))
for(i in 1: nrow(sites)){
  hq_pred[i,] = sm_forward_evap(sites$map.wc[i], sites$mat.wc[i], sites$hqp.frac[i], sites$hqt.offset[i], sites$Ra[i], 280)
  
}
hq_pred$Site = sites$Site

dq_pred = data.frame(d13C=numeric(0), d18O=numeric(0), d13C_sd=numeric(0), d18O_sd=numeric(0))
for(i in 1: nrow(sites)){
  dq_pred[i,] = sm_forward_evap(sites$map.wc[i], sites$mat.wc[i], sites$dqp.frac[i], sites$dqt.offset[i], sites$Ra[i], 280)
  
}
dq_pred$Site = sites$Site


## Plots w/ MAP colors
hq.comp = merge.data.frame(hq_pred, sites, by.x = "Site", by.y = "Site", all.x=TRUE)
dq.comp = merge.data.frame(dq_pred, sites, by.x = "Site", by.y = "Site", all.x=TRUE)

hq.comp = merge.data.frame(hq.comp, data.aves, by.x = "Site", by.y = "Group.1", all.x=TRUE)
dq.comp = merge.data.frame(dq.comp, data.aves, by.x = "Site", by.y = "Group.1", all.x=TRUE)

View(hq.comp)

library(RColorBrewer)

c = ceiling(((hq.comp$map.wc - 100) / 700) * 7)

pal_map = brewer.pal(7, "Blues")

c_mat = ceiling(((hq.comp$mat.wc + 5) / 20 * 5))

pal_mat = brewer.pal(6, "YlOrRd")

jpeg("Pre_Validation.jpg", units="in", res=300, width=8.5, height=7.5)

layout(matrix(c(1,2,3,4), 2, 2, byrow=T), heights=c(1,1.1,1,1.1), widths=c(1,1.05,1,1.05))
par(mar=c(2,5,1,0))
plot(hq.comp$d13C.measured, hq.comp$d13C, pch=16, col=pal_map[c], xlim=c(-12,2.1), ylim=c(-12,2.1), main="", cex = 1.25,
     xlab="", cex.axis=1.25, cex.lab=1.25,
     ylab=expression(paste("Predicted ",delta^{13}, "C (\u2030)")),
     xaxt='n')
abline(0,1)
points(hq.comp$d13C.measured, hq.comp$d13C, pch=1)
text(-12,2,"a", cex=1.5)

par(mar=c(2,5,1,1))
plot(hq.comp$d18O.measured, hq.comp$d18O, pch=16, col=pal_map[c], xlim=c(-16,0), ylim=c(-16,0),main="",cex = 1.25,
     xlab="",  cex.axis=1.25, cex.lab=1.25,
     ylab=expression(paste("Predicted ",delta^{18}, "O (\u2030)")),
     xaxt='n')
abline(0,1)
points(hq.comp$d18O.measured, hq.comp$d18O, pch=1)
text(-16,0,"b", cex=1.5)

par(mar=c(4,5,0,0))
plot(dq.comp$d13C.measured, dq.comp$d13C, pch=16, col=pal_map[c], xlim=c(-12,2.1), ylim=c(-12,2.1), main="",cex = 1.25,
     xlab=expression(paste("Observed ",delta^{13}, "C (\u2030)")),  cex.axis=1.25, cex.lab=1.25,
     ylab=expression(paste("Predicted ",delta^{13}, "C (\u2030)")))
abline(0,1)
points(dq.comp$d13C.measured, dq.comp$d13C, pch=1)
text(-12,2,"c", cex=1.5)

par(mar=c(4,5,0,1))
plot(dq.comp$d18O.measured, dq.comp$d18O, pch=16, col=pal_map[c], xlim=c(-16,0), ylim=c(-16,0), main="",cex = 1.25,
     xlab=expression(paste("Observed ",delta^{18}, "O (\u2030)")),  cex.axis=1.25, cex.lab=1.25,
     ylab=expression(paste("Predicted ",delta^{18}, "O (\u2030)")))
abline(0,1)
points(dq.comp$d18O.measured, dq.comp$d18O, pch=1)
text(-16,0,"d", cex=1.5)

legend("bottomright", title = expression(paste("P"[a]," (mm)")), cex=1.05, fill=pal_map, legend=c("100 - 200", "200 - 300", "300 - 400", "400 - 500", "500 - 600", "600 - 700", "700 - 800"))

dev.off()

## Basic stats
pre_lm_hot_C <- lm()







