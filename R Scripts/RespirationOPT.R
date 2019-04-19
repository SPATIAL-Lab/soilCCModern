sm_optimizer_r = function(MAP, MAT, P_seas, T_seas, pCO2, Ra, rr){
  
  deltaA = rnorm(nsynth, -6.5, 0.3)
  pores = rnorm(nsynth, 0.4, 0.1)
  tort = rnorm(nsynth, 0.6, 0.1)
  
  #Solar radiation, here fixed
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
  CQT_K <- CQT + 273.15
  
  #Soil Temps
  d <- sqrt((2*0.0007)/((2*3.1415/3.154e7)*0.3))
  t <- ifelse(CQT > MAT + 3, 1, ifelse(CQT < MAT - 3, 3, 4))
  T_soil <- MAT + (hqt.offset * sin((2*3.1415/4) * t - z/d) / exp(z/d)) 
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
  
  #Depth to carbonate, now gamma dist
  z_mean <- MAP*0.0925 + 13.4  #top of Bk equation based on Retallack 2005 data
  z_thick = abs(CMP_mm - MAP/12) * 0.74 + 17.3 #thickness of Bk, CMP as proxy for seasonality 
  z_mean = z_mean + z_thick/2 #find middle of Bk in cm
  theta = 20^2/z_mean  #gamma scale parameter, using 20cm as variance from Retallack fit
  k = z_mean / theta #gamma shape parameter
  z = rgamma(nsynth, shape = k, scale = theta)
  z <- pmin(z, 100) 
  z_m <- z * 0.01
  
  #Respiration rate, now gamma dist
  R_month_m <- 1.25 * exp(0.055 * CQT) * CMP_cm / (4.78 + CMP_cm)  #Raich 2002, gC/m2day
  R_month_m <- R_month_m * rr
  theta = (R_month_m*0.5)^2/R_month_m #gamma scale parameter, using mean residual of 50% based on Raich validation data 
  k = R_month_m / theta #gamma shape parameter
  R_month = rgamma(nsynth, shape = k, scale = theta) #lets use gamma for these quants bounded at zero....
  R_month = R_month / (12.01 * 100^2)  #molC/cm2month
  R_sec <- R_month / (24 * 3600)  #molC/cm2s
  R_sec = R_sec / L / pores #molC/cm3s
  
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
  
  #CO2 Diffusion coefficient
  DIFC = FAP * tort * 0.1369 * (T_soil_K / 273.15) ^ 1.958
  
  #Water limitation of discriminaton, Diefendorf
  W_m <- 22.65 - (1.2 * (MAP + 975)) / (27.2 + 0.04 * (MAP + 975))
  W = rnorm(nsynth, W_m, 0.5)
  
  #CO2 effect on discrimination, Schubert
  deltaP_pCO2_m <- 28.26 * 0.35 * (pCO2 + 15) / (28.26 + 0.35 * (pCO2 + 15))
  deltaP_pCO2 = rnorm(nsynth, deltaP_pCO2_m, 0.5)
  
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
  
  dC_Carb <- (R_Carb / RC.vpdb - 1) * 1000
  dO_Carb <- (R_O_Carb / RO.vpdb - 1) * 1000
  
  return(median(dC_Carb))
}

rr = seq(0.01, 1, 0.01)

parms_hq = data.frame(rr = rr, rmse = numeric(100))

## HQ

for(j in 1:nrow(parms_hq)){
  opt = numeric()
  for(i in 1: nrow(sites)){
    opt[i] = sm_optimizer_r(sites$map.wc[i], sites$mat.wc[i], sites$hqp.frac[i], sites$hqt.offset[i], 280, sites$Ra[i], parms_hq$rr[j])
  }
  
  opt = data.frame(Site = sites$Site, d13C = opt)
  opt = merge.data.frame(opt, data.comp, by.x = "Site", by.y = "Site", all.x=TRUE)
  
  mse = (opt$d13C - opt$d13C.measured)^2
  mse = mean(mse)
  parms_hq$rmse[j] = sqrt(mse)
}

rr_opt_hq = parms_hq[which.min(parms_hq$rmse),"rr"]

## DQ
  
parms_dq = data.frame(rr = rr, rmse = numeric(100))

for(j in 1:nrow(parms_dq)){
  opt = numeric()
  for(i in 1: nrow(sites)){
    opt[i] = sm_optimizer_r(sites$map.wc[i], sites$mat.wc[i], sites$dqp.frac[i], sites$dqt.offset[i], 280, sites$Ra[i], parms_dq$rr[j])
  }
  
  opt = data.frame(Site = sites$Site, d13C = opt)
  opt = merge.data.frame(opt, data.comp, by.x = "Site", by.y = "Site", all.x=TRUE)
  
  mse = (opt$d13C - opt$d13C.measured)^2
  mse = mean(mse)
  parms_dq$rmse[j] = sqrt(mse)
}

rr_opt_dq = parms_dq[which.min(parms_dq$rmse),"rr"]


jpeg("Resp_OPT.jpg", units="in", res=300, width=10, height=4.68)

layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths = c(1.25,1))

par(mar=c(4,5,1,1))
plot(parms_hq$rmse ~ parms_hq$rr, cex = 1.5, cex.axis=1.3, cex.lab=1.3, ylim=c(1,10.5), type="n", main = "", ylab = expression(paste("RMSE (\u2030)")), xlab = "Fraction of Estimated Respiration")
lines(parms_hq$rmse ~ parms_hq$rr)
text(1,10.4,"a", cex=1.5)
par(mar=c(4,0,1,1))
plot(parms_dq$rmse ~ parms_dq$rr, cex = 1.5, cex.axis=1.3, cex.lab=1.3, type="n", yaxt='n', ylim=c(1,10.5), main = "", ylab = "RMSE", xlab = "Fraction of Estimated Respiration")
lines(parms_dq$rmse ~ parms_dq$rr)
text(1,10.4,"b", cex=1.5)

dev.off()



