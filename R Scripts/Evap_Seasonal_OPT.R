############### Forward model function for use in sensitivity testing
sm_optimizer_evap = function(MAP, MAT, P_seas, T_seas, pCO2, Ra, spre, esw){

  deltaA = rnorm(nsynth, -6.5, 0.3)
  pores = rnorm(nsynth, 0.4, 0.1)
  tort = rnorm(nsynth, 0.7, 0.1)
  
  #Solar radiation, here based on annual avg at latitude of sites
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
  T_OOS <- (MAT * 4 -  CQT) / 3
  
  #Soil Temps
  d <- sqrt((2*0.0007)/((2*3.1415/3.154e7)*0.3))
  t <- ifelse(CQT > MAT + 3, 0.3, ifelse(CQT < MAT - 3, 0.8, 0.05))
  T_soil <- MAT + (hqt.offset * sin((2*3.1415) * t - z/d) / exp(z/d)) 
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
  
  #Precipitation O isotope ratios 
  dO_P_OOS <- -14.76 + 0.56 * T_OOS  #Relevant precip is spre % from CQ
  dO_P_PCQ <- -14.76 + 0.56 * CQT
  dO_P_m <- (dO_P_PCQ * P_seas + dO_P_OOS * (1 - spre) * (1 - P_seas)) / (P_seas + (1 - spre) * (1 - P_seas))
  dO_P = rnorm(nsynth, dO_P_m, 1.46)
  R_O_P = (dO_P / 1000 + 1) * RO.vsmow
  
  #Atmospheric vapor O isotope ratios
  A_atmP <- 2.71828 ^ ((5.9702e6 / CQT_K ^ 2 - 3.2801e4 / CQT_K + 52.227) / 1000)
  R_O_atm <- R_O_P / A_atmP
  dO_atm_m <- (R_O_atm / RO.vsmow - 1) * 1000
  dO_atm = rnorm(nsynth, dO_atm_m, 1)
  
  #Depth to carbonate, now gamma dist
  z_mean <- MAP*0.0925 + 13.4  #top of Bk equation based on Retallack 2005 data
  z_thick = abs(CMP_mm - MAP/12) * 0.74 + 17.3 #thickness of Bk, CMP as proxy for seasonality 
  z_mean = z_mean + z_thick/2 #find middle of Bk in cm
  theta = 20^2/z_mean  #gamma scale parameter, using 20cm as variance from Retallack fit
  k = z_mean / theta #gamma shape parameter
  z = rgamma(nsynth, shape = k, scale = theta)
  z <- pmin(z, 100) 
  z_m <- z * 0.01
  
  #Potential ET
  ETP_D_m <- ifelse (RH < 50, 0.013 * (CQT / (CQT + 15)) * (23.8856 * Rs + 50) * (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50))
  ETP_M <- ETP_D * 30  #mm/month
  
  #Actual ET
  ETA_var = rnorm(nsynth, 1, 0.2) #This noise parmeter limits ETA<CMP_mm but allows variation around ETP, as observed
  ETA = CMP_mm*3 * (1 / (sqrt(1 + (1 / ((ETP_M / (CMP_mm*3)) * ETA_var)) ^ 2))) #AET in mm/month from Budyko curve
  #here scaled eta to quarter precip, assuming potential carry-over
  
  #Free air porosity
  #Have updated, now scales volumetrically w/ excess precipitation relative to pore space
  CQP_soil <- (MAP * (P_seas + (1 - spre) * (1 - P_seas))) / 3
  FAP <- pmin((pores - (CQP_soil - ETA)/(L * 10 * pores)), pores-0.05)
  FAP = pmax(FAP,0.01) #dimensionless
  
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
  DIFO <- 1.637e-8 * (T_soil_K / 216.25 - 1) ^ 2.074 * (pores-FAP) * tort   ## soil water fraction
  z_i <- DIFO / E_s #mean penetration depth of evap, in m
  
  #Soil water O isotopes
  DRF <- 1 + 0.8 * (1 / 0.9723 - 1)
  R_O_surface <- ((1 - h) * DRF * R_O_P + h * R_O_atm) / (1 / A_atmP)
  R_O_soil_evap <- ((R_O_surface - R_O_P) * 2.71828 ^ (-z_m / z_i)) + R_O_P
  R_O_soil <- R_O_soil_evap * esw + R_O_P * (1 - esw)  #soil water is esw % evaporated fraction
  R_O_soil <- ifelse(is.na(R_O_soil), R_O_P, R_O_soil)
  #Soil carbonate O isotopes
  A_O <- 2.71828 ^ ((2.78e6 / T_soil_K ^ 2 - 2.89) / 1000)
  R_O_Carb <- R_O_soil * A_O
  dO_Carb <- (R_O_Carb / RO.vpdb - 1) * 1000
  
  dat = c(median(dO_Carb))
  return(dat)
}

s = seq(0,0.95,0.05)
ss = seq(0, 1-0.05/20, 0.05/20)
ss = ss*20
ss = trunc(ss)
ss = ss/20

parms_hq = data.frame(spres = ss, esws = rep(s, 20), rmse = numeric(400), d18O = numeric(400), R_O_soil = numeric(400))

for(j in 1:nrow(parms_hq)){
  opt = numeric()
  for(i in 1: nrow(sites)){
    opt[i] = sm_optimizer_evap(sites$map.wc[i], sites$mat.wc[i], sites$hqp.frac[i], sites$hqt.offset[i], 280, sites$Ra[i], parms$spres[j], parms$esws[j])
  }
  
  opt = data.frame(Site = sites$Site, d18O = opt)
  opt = merge.data.frame(opt, data.comp, by.x = "Site", by.y = "Site", all.x=TRUE)
  
  mse = (opt$d18O - opt$d18O.measured)^2
  mse = mean(mse)
  parms_hq$rmse[j] = sqrt(mse)
  parms_hq$d18O[j] = opt$d18O
  parms_hq$R_O_Soil[j] = opt$R_O_P
}

min(parms_hq$rmse)
View(parms_hq)
rmses_hq = matrix(parms_hq$rmse, 20, 20)
rmses_hq = rmses_hq[c(20:1),]
rmses.rast_hq = raster(rmses_hq, xmn=0, xmx=0.95, ymn=0, ymx=0.95)


parms_dq = data.frame(spres = ss, esws = rep(s, 20), rmse = numeric(400))

for(j in 1:nrow(parms_dq)){
  opt = numeric()
  for(i in 1: nrow(sites)){
    opt[i] = sm_optimizer_evap(sites$map.wc[i], sites$mat.wc[i], sites$hqp.frac[i], sites$hqt.offset[i], 280, sites$Ra[i], parms$spres[j], parms$esws[j])
  }
  
  opt = data.frame(Site = sites$Site, d18O = opt)
  opt = merge.data.frame(opt, data.comp, by.x = "Site", by.y = "Site", all.x=TRUE)
  
  mse = (opt$d18O - opt$d18O.measured)^2
  mse = mean(mse)
  parms_dq$rmse[j] = sqrt(mse)
}

View(parms_dq)
min(parms_dq$rmse)
rmses_dq = matrix(parms_dq$rmse, 20, 20)
rmses_dq = rmses_dq[c(20:1),]
rmses.rast_dq = raster(rmses_dq, xmn=0, xmx=0.95, ymn=0, ymx=0.95)

jpeg("Evap_OPT.jpg", units="in", res=300, width=10, height=4.68)
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(1.001,1))

par(mar=c(5.5, 5, 0, 0))
plot(rmses.rast_hq, main="", cex.lab=1.3, cex.axis=1.3, xlab="% Seasonal Rainfall", ylab="% Evaporated Water", xlim=c(0,1), ylim=c(0,1), zlim=c(2.5,6.5), legend = F)

par(mar=c(5.5, 0, 0, 7))
plot(rmses.rast_dq, yaxt='n', cex.lab=1.3, cex.axis=1.3, main="", xlab="% Seasonal Rainfall", ylab="", xlim=c(0,1), ylim=c(0,1), zlim=c(2.5,6.5))  #now need to make a nice plot...
mtext("RMSE, (\u2030)", 4, line=0.5, cex=1.3)

dev.off()



