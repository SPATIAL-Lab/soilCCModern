sm_forward_evap = function(MAP, MAT, P_seas, T_seas, Ra, pCO2){
  
  deltaA = rnorm(nsynth, -6.5, 0.3)
  pores = rnorm(nsynth, 0.4, 0.1)
  tort = rnorm(nsynth, 0.7, 0.1)
  esw = 1
  
  #Solar radiation, here based on yearly average at latitude
  Rs = Ra * 0.16 * sqrt(12)
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
  
  #Depth to carbonate, now gamma dist
  z_mean <- MAP*0.0925 + 13.4  #top of Bk equation based on Retallack 2005 data
  z_thick = abs(CMP_mm - MAP/12) * 0.74 + 17.3 #thickness of Bk, CMP as proxy for seasonality 
  z_mean = z_mean + z_thick/2 #find middle of Bk in cm
  theta = 20^2/z_mean  #gamma scale parameter, using 20cm as variance from Retallack fit
  k = z_mean / theta #gamma shape parameter
  z = rgamma(nsynth, shape = k, scale = theta)
  z <- pmin(z, 100)
  
  #Soil Temperature
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
  
  # Precipitation O isotope ratios 
  dO_P_m <- -14.76 + 0.56 * MAT
  dO_P = rnorm(nsynth, dO_P_m, 1.46)
  R_O_P = (dO_P / 1000 + 1) * RO.vsmow
  
 
  
  #Respiration rate, now gamma dist
  R_day_m <- 1.25 * exp(0.0545 * CQT) * CMP_cm / (4.259 + CMP_cm)  #Raich 2002, gC/m2day
  theta = (R_day_m*0.5)^2/R_day_m #gamma scale parameter, using mean residual of 50% based on Raich validation data 
  k = R_day_m / theta #gamma shape parameter
  R_day = rgamma(nsynth, shape = k, scale = theta) #lets use gamma for these quants bounded at zero....
  R_day = R_day / (12.01 * 100^2)  #molC/cm2day
  R_sec <- R_day / (24 * 3600)  #molC/cm2s
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
  
  #Soil carbonate C isotopes

  A_CO2_Carb <- 2.71828 ^ (-2.988e3 / T_soil_K ^ 2 + 7.6663 / T_soil_K - 0.0024612)
  R_Soil <- (dC_Soil / 1000 + 1) * RC.vpdb
  R_Carb <- R_Soil / A_CO2_Carb
  
  #Soil water evaporation, now beta dist
  e_mean = 0.06  #evap is 6% of total ET
  e_var = 0.04^2
  size = e_mean*(1-e_mean)/e_var - 1
  alpha = e_mean * size
  beta = (1-e_mean) * size
  E = rbeta(nsynth, alpha, beta) * ETA #mm/month
  E = pmax(E, 1) #mm/month
  
  #Soil water diffusion evaporation balance
  E_s <- E / (1000 * 30 * 24 * 3600) #evaporation in m/sec
  DIFO <- 1.637e-8 * (T_soil_K / 216.25 - 1) ^ 2.074 * (pores - FAP) * tort   ## should be soil water fraction, 
  ## pores - FAP. units: m2/sec. However, the the paper assumes total saturation, where FAP = 0
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

## Try idea: Low MAP is dry quarter, high MAP is hot quarter. 
cut = seq(110, 700, 2)
cutoff = data.frame(cut, rmse_c = numeric(296), rmse_o = numeric(296))

for(j in 1:nrow(cutoff)){
  sites_dry <- subset.data.frame(sites, map.wc <= cutoff$cut[j])
  sites_wet <- subset.data.frame(sites, map.wc > cutoff$cut[j])
  
  hq_pred_wet = data.frame(d13C=numeric(0), d18O=numeric(0), d13C_sd=numeric(0), d18O_sd=numeric(0))
  for(i in 1: nrow(sites_wet)){
    hq_pred_wet[i,] = sm_forward_evap(sites_wet$map.wc[i], sites_wet$mat.wc[i], sites_wet$hqp.frac[i], sites_wet$hqt.offset[i], sites_wet$Ra[i], 280)
    
  }
  hq_pred_wet$Site = sites_wet$Site
  
  opt_c_wet = data.frame(Site = sites_wet$Site, d13C = hq_pred_wet$d13C)
  opt_o_wet = data.frame(Site = sites_wet$Site, d18O = hq_pred_wet$d18O)
  opt_c_wet = merge.data.frame(opt_c_wet, data.comp, by.x = "Site", by.y = "Site", all.x=TRUE)
  opt_o_wet = merge.data.frame(opt_o_wet, data.comp, by.x = "Site", by.y = "Site", all.x=TRUE)
  
  dq_pred_dry = data.frame(d13C=numeric(0), d18O=numeric(0), d13C_sd=numeric(0), d18O_sd=numeric(0))
  for(i in 1: nrow(sites_dry)){
    dq_pred_dry[i,] = sm_forward_evap(sites_dry$map.wc[i], sites_dry$mat.wc[i], sites_dry$dqp.frac[i], sites_dry$dqt.offset[i], sites_dry$Ra[i], 280)
    
  }
  
  dq_pred_dry$Site = sites_dry$Site
  
  opt_c_dry = data.frame(Site = sites_dry$Site, d13C = dq_pred_dry$d13C)
  opt_o_dry = data.frame(Site = sites_dry$Site, d18O = dq_pred_dry$d18O)
  opt_c_dry = merge.data.frame(opt_c_dry, data.comp, by.x = "Site", by.y = "Site", all.x=TRUE)
  opt_o_dry = merge.data.frame(opt_o_dry, data.comp, by.x = "Site", by.y = "Site", all.x=TRUE)
  
  opt_c <- rbind(opt_c_wet, opt_c_dry)
  opt_o <- rbind(opt_o_wet, opt_o_dry)
  
  mse_c = (opt_c$d13C - opt_c$d13C.measured)^2
  mse_o = (opt_o$d18O - opt_o$d18O.measured)^2
  mse_c = mean(mse_c)
  mse_o = mean(mse_o)
  
  cutoff$rmse_c[j] <- sqrt(mse_c)
  cutoff$rmse_o[j] <- sqrt(mse_o)
  
}

hq_pred_wet = data.frame(d13C=numeric(0), d18O=numeric(0), d13C_sd=numeric(0), d18O_sd=numeric(0))
for(i in 1: nrow(sites_wet)){
  hq_pred_wet[i,] = sm_forward_evap(sites_wet$map.wc[i], sites_wet$mat.wc[i], sites_wet$hqp.frac[i], sites_wet$hqt.offset[i], sites_wet$Ra[i], 280)
  
}
hq_pred_wet$Site = sites_wet$Site

dq_pred_dry = data.frame(d13C=numeric(0), d18O=numeric(0), d13C_sd=numeric(0), d18O_sd=numeric(0))
for(i in 1: nrow(sites_dry)){
  dq_pred_dry[i,] = sm_forward_evap(sites_dry$map.wc[i], sites_dry$mat.wc[i], sites_dry$dqp.frac[i], sites_dry$dqt.offset[i], sites_dry$Ra[i], 280)
  
}
dq_pred_dry$Site = sites_dry$Site

pred <- rbind(dq_pred_dry, hq_pred_wet)
pred <- merge.data.frame(pred, sites, by.x = "Site")
pred <- merge.data.frame(pred, data.aves, by.x = "Site", by.y = "Group.1", all.x=TRUE)



cutoff$rmse <- (cutoff$rmse_c + cutoff$rmse_o) / 2

jpeg("Cutoff.jpg", res=300, units="in", height=8, width=8)
par(mar=c(5,5,1,1))
plot(cutoff$rmse, cutoff$cut, type = "l", xlab = "Mean RMSE", ylab = expression(paste("Cutoff P"[a]," (mm)")), cex.lab = 1.25, cex.axis= 1.25)
dev.off()

plot(cutoff$rmse_c, cutoff$cut, type="l")
plot(cutoff$rmse_o, cutoff$cut, type="l")

## Minimum RMSE (combined) has a 280 mm cutoff


c = ceiling(((pred$map.wc - 100) / 700) * 7)

pal_map = brewer.pal(7, "Blues")

c_mat = ceiling(((pred$mat.wc + 5) / 20 * 5))

pal_mat = brewer.pal(6, "YlOrRd")

layout(matrix(c(1,2), 2, 2, byrow=T))

plot(pred$d13C.measured, pred$d13C, pch=16, col=pal_map[c], xlim=c(-12,1.5), ylim=c(-12,1.5), main="", cex = 1.25,
     xlab=expression(paste("Observed ",delta^{13}, "C (\u2030)")),
     ylab=expression(paste("Predicted ",delta^{13}, "C (\u2030)")))
abline(0,1)
points(pred$d13C.measured, pred$d13C, pch=1)

plot(pred$d18O.measured, pred$d18O, pch=16, col=pal_mat[c_mat], xlim=c(-17,4), ylim=c(-17,4),main="",cex = 1.25,
     xlab=expression(paste("Observed ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Predicted ",delta^{18}, "O (\u2030)")))
abline(0,1)
points(pred$d18O.measured, pred$d18O, pch=1)
