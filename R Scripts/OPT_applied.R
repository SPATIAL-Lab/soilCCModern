sm_forward_opt_hq = function(MAP, MAT, P_seas, T_seas, Ra, pCO2){
  
  deltaA = rnorm(nsynth, -6.5, 0.3)
  pores = rnorm(nsynth, 0.4, 0.1)
  tort = rnorm(nsynth, 0.7, 0.1)
  
  #Solar radiation, by latitude
  Rs = Ra * sqrt(12) * 0.16
  DIF.ratio = 1.004443
  
  #Isotope ratio constants
  RC.vpdb = 0.011237
  RO.vsmow = 0.0020052
  RO.vpdb = 0.002067181
  
  #Basal depth
  L = 100
  
  #Climate stuff
  CQP <- MAP * P_seas
  CQP = max(CQP, 3)
  CMP_mm <- CQP / 3
  CMP_cm <- CQP / 30
  CQT <- MAT + T_seas
  CQT_K <- CQT + 273
  
  #Soil temps
  t <- 1
  d <- sqrt((2*0.0007)/((2*3.1415/3.154e7)*0.3))
  T_soil <- MAT + (T_seas * sin((2*3.1415/4) * t - z/d) / exp(z/d)) 
  T_soil_K <- T_soil + 273.15
  
  #Convert pCO2 to units mol/cm^3
  pCO2_mcc = pCO2 / (0.08206 * T_soil_K * 10^9)  #mol/cm^3
  
  #Relative humidity, now beta distribution
  h_m <- 0.25 + 0.7 * (CQP / 900)
  h_var = 0.05^2
  size = h_m*(1-h_m)/h_var - 1
  alpha = h_m * size
  beta = (1-h_m) * size
  h = rbeta(nsynth, alpha, beta)
  RH <- h * 100
  
  # Precipitation O isotope ratios 
  dO_P_m <- -14.76 + 0.56 * MAT
  dO_P = rnorm(nsynth, dO_P_m, 1.46)
  R_O_P = (dO_P / 1000 + 1) * RO.vsmow
  
  #Depth to carbonate, now gamma dist
  z_mean <- MAP*0.0925 + 13.4  #top of Bk equation based on Retallack 2005 data
  z_thick = abs(CMP_mm - MAP/12) * 0.74 + 17.3 #thickness of Bk, CMP as proxy for seasonality 
  z_mean = z_mean + z_thick/2 #find middle of Bk in cm
  theta = 20^2/z_mean  #gamma scale parameter, using 20cm as variance from Retallack fit
  k = z_mean / theta #gamma shape parameter
  z = rgamma(nsynth, shape = k, scale = theta)
  z <- pmin(z, 100) 
  
  #Respiration rate, now gamma dist
  R_day_m <- 1.25 * exp(0.0545 * CQT) * CMP_cm / (4.259 + CMP_cm)  #Raich 2002, gC/m2day
  R_day_m <- R_day_m * rr_opt_hq # Hot quarter respiration rate optimization - 39% of estimated
  theta = (R_day_m*0.5)^2/R_day_m #gamma scale parameter, using mean residual of 50% based on Raich validation data 
  k = R_day_m / theta #gamma shape parameter
  R_day = rgamma(nsynth, shape = k, scale = theta) #lets use gamma for these quants bounded at zero....
  R_day = R_day / (12.01 * 100^2)  #molC/cm2day
  R_sec <- R_day / (24 * 3600)  #molC/cm2s
  R_sec = R_sec / L / pores #molC/cm3s
  
  #
  
  #Potential ET
  ETP_D_m <- ifelse (RH < 50, 0.013 * (CQT / (CQT + 15)) * (23.8856 * Rs + 50) * (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50))
  ETP_D = rnorm(nsynth, ETP_D_m, 0.2)  #PET in mm/day, Turc 1961
  ETP_M <- ETP_D * 30  #mm/month
  
  #Actual ET
  ETA_var = rnorm(nsynth, 1, 0.2) #This noise parmeter limits ETA<CMP_mm but allows variation around ETP, as observed
  ETA = CMP_mm*3 * (1 / (sqrt(1 + (1 / ((ETP_M / (CMP_mm*3)) * ETA_var)) ^ 2))) #AET in mm/month from Budyko curve
  #here scaled eta to quarter precip, assuming potential carry-over
  
  #Free air porosity
  #Have updated, now scales volumetrically w/ excess precipitation relative to pore space
  FAP <- pmin((pores - (CMP_mm - ETA)/(L*10*pores)), pores - 0.05)
  FAP = pmax(FAP,0.01) #dimensionless
  
  #CO2 Diffusion coefficients
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
  

  A_CO2_Carb <- 2.71828 ^ (-2.988e3 / T_soil_K ^ 2 + 7.6663 / T_soil_K - 0.0024612)
  R_Soil <- (dC_Soil / 1000 + 1) * RC.vpdb
  R_Carb <- R_Soil / A_CO2_Carb
  
  #Soil carbonate O isotopes
  A_O <- 2.71828 ^ ((2.78e6 / T_soil_K ^ 2 - 2.89) / 1000)
  R_O_Carb <- R_O_P * A_O
  
  dC_Carb <- (R_Carb / RC.vpdb - 1) * 1000 
  dO_Carb <- (R_O_Carb / RO.vpdb - 1) * 1000
  
  dat = c(median(dC_Carb), median(dO_Carb), sd(dC_Carb), sd(dO_Carb))
  
  return(dat)
}

sm_forward_opt_dq = function(MAP, MAT, P_seas, T_seas, Ra, pCO2){
  
  deltaA = rnorm(nsynth, -6.5, 0.3)
  pores = rnorm(nsynth, 0.4, 0.1)
  tort = rnorm(nsynth, 0.7, 0.1)
  
  #Solar radiation, by latitude
  Rs = Ra * sqrt(12) * 0.16
  DIF.ratio = 1.004443
  
  #Isotope ratio constants
  RC.vpdb = 0.011237
  RO.vsmow = 0.0020052
  RO.vpdb = 0.002067181
  
  #Basal depth
  L = 100
  
  #Climate stuff
  CQP <- MAP * P_seas
  CQP = max(CQP, 3)
  CMP_mm <- CQP / 3
  CMP_cm <- CQP / 30
  CQT <- MAT + T_seas
  CQT_K <- CQT + 273
  
  #Convert pCO2 to units mol/cm^3
  pCO2_mcc = pCO2 / (0.08206 * CQT_K * 10^9)  #mol/cm^3
  
  #Relative humidity, now beta distribution
  h_m <- 0.25 + 0.7 * (CQP / 900)
  h_var = 0.05^2
  size = h_m*(1-h_m)/h_var - 1
  alpha = h_m * size
  beta = (1-h_m) * size
  h = rbeta(nsynth, alpha, beta)
  RH <- h * 100
  
  # Precipitation O isotope ratios 
  dO_P_m <- -14.76 + 0.56 * MAT
  dO_P = rnorm(nsynth, dO_P_m, 1.46)
  R_O_P = (dO_P / 1000 + 1) * RO.vsmow
  
  #Depth to carbonate, now gamma dist
  z_mean <- MAP*0.0925 + 13.4  #top of Bk equation based on Retallack 2005 data
  z_thick = abs(CMP_mm - MAP/12) * 0.74 + 17.3 #thickness of Bk, CMP as proxy for seasonality 
  z_mean = z_mean + z_thick/2 #find middle of Bk in cm
  theta = 20^2/z_mean  #gamma scale parameter, using 20cm as variance from Retallack fit
  k = z_mean / theta #gamma shape parameter
  z = rgamma(nsynth, shape = k, scale = theta)
  z <- pmin(z, 100) 
  
  #Respiration rate, now gamma dist
  R_day_m <- 1.25 * exp(0.0545 * CQT) * CMP_cm / (4.259 + CMP_cm)  #Raich 2002, gC/m2day
  R_day_m <- R_day_m * rr_opt_dq # Dry site respiration rate optimization - 46% estimated
  theta = (R_day_m*0.5)^2/R_day_m #gamma scale parameter, using mean residual of 50% based on Raich validation data 
  k = R_day_m / theta #gamma shape parameter
  R_day = rgamma(nsynth, shape = k, scale = theta) #lets use gamma for these quants bounded at zero....
  R_day = R_day / (12.01 * 100^2)  #molC/cm2day
  R_sec <- R_day / (24 * 3600)  #molC/cm2s
  R_sec = R_sec / L / pores #molC/cm3s
  
  #
  
  #Potential ET
  ETP_D_m <- ifelse (RH < 50, 0.013 * (CQT / (CQT + 15)) * (23.8856 * Rs + 50) * (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50))
  ETP_M <- ETP_D * 30  #mm/month
  
  #Actual ET
  ETA_var = rnorm(nsynth, 1, 0.2) #This noise parmeter limits ETA<CMP_mm but allows variation around ETP, as observed
  ETA = CMP_mm*3 * (1 / (sqrt(1 + (1 / ((ETP_M / (CMP_mm*3)) * ETA_var)) ^ 2))) #AET in mm/month from Budyko curve
  #here scaled eta to quarter precip, assuming potential carry-over
  
  #Free air porosity
  #Have updated, now scales volumetrically w/ excess precipitation relative to pore space
  FAP <- pmin((pores - (CMP_mm - ETA)/(L*10*pores)), pores-0.05)
  FAP = pmax(FAP,0.01) #dimensionless
  
  #CO2 Diffusion coefficients
  DIFC = FAP * tort * 0.1369 * (CQT_K / 273.15) ^ 1.958
  
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
  
  #Soil carbonate C isotopes
  t <- 4
  d <- sqrt((2*0.0007)/((2*3.1415/3.154e7)*0.3))
  T_soil <- MAT + (hqt.offset * sin((2*3.1415/4) * t - z/d) / exp(z/d)) 
  T_soil_K <- T_soil + 273.15
  A_CO2_Carb <- 2.71828 ^ (-2.988e3 / T_soil_K ^ 2 + 7.6663 / T_soil_K - 0.0024612)
  R_Soil <- (dC_Soil / 1000 + 1) * RC.vpdb
  R_Carb <- R_Soil / A_CO2_Carb
  
  #Soil carbonate O isotopes
  A_O <- 2.71828 ^ ((2.78e6 / T_soil_K ^ 2 - 2.89) / 1000)
  R_O_Carb <- R_O_P * A_O
  
  dC_Carb <- (R_Carb / RC.vpdb - 1) * 1000 
  dO_Carb <- (R_O_Carb / RO.vpdb - 1) * 1000
  
  dat = c(median(dC_Carb), median(dO_Carb), sd(dC_Carb), sd(dO_Carb))
  
  return(dat)
}



## Post- validation after optimization

hq_pred = data.frame(d13C=numeric(0), d18O=numeric(0), d13C_sd=numeric(0), d18O_sd=numeric(0))
for(i in 1: nrow(sites)){
  hq_pred[i,] = sm_forward_opt_hq(sites$map.wc[i], sites$mat.wc[i], sites$hqp.frac[i], sites$hqt.offset[i], sites$Ra[i], 280)
  
}
hq_pred$Site = sites$Site

dq_pred = data.frame(d13C=numeric(0), d18O=numeric(0), d13C_sd=numeric(0), d18O_sd=numeric(0))
for(i in 1: nrow(sites)){
  dq_pred[i,] = sm_forward_opt_dq(sites$map.wc[i], sites$mat.wc[i], sites$dqp.frac[i], sites$dqt.offset[i], sites$Ra[i], 280)
  
}
dq_pred$Site = sites$Site

hq.comp = merge.data.frame(hq_pred, sites, by.x = "Site", by.y = "Site", all.x=TRUE)
dq.comp = merge.data.frame(dq_pred, sites, by.x = "Site", by.y = "Site", all.x=TRUE)

hq.comp = merge.data.frame(hq.comp, data.aves, by.x = "Site", by.y = "Group.1", all.x=TRUE)
dq.comp = merge.data.frame(dq.comp, data.aves, by.x = "Site", by.y = "Group.1", all.x=TRUE)



#add it to the predictions and plot
library(RColorBrewer)

c = ceiling(((hq.comp$map.wc - 100) / 700) * 7)

pal_map = brewer.pal(7, "Blues")

c_mat = ceiling(((hq.comp$mat.wc + 5) / 20 * 5))

pal_mat = brewer.pal(6, "YlOrRd")

jpeg("Post_Validation.jpg", units="in", res=300, width=8.5, height=7.5)

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

legend("bottomright", title = "MAP (mm)", cex=1.05, fill=pal_map, legend=c("100 - 200", "200 - 300", "300 - 400", "400 - 500", "500 - 600", "600 - 700", "700 - 800"))

dev.off()

View(hq.comp)
View(dq.comp)

