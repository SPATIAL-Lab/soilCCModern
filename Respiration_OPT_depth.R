sm_optimizer_r.depth = function(Pa, Tma, PfPCQ, TmPCQ_min_a, Ra, TmWQ_min_a, dryq, rr, z){
  
  #Pre-industrial pCO2 and carbon isotope value of atmospheric CO2
  d13C_atmCO2 = rnorm(nsynth, -6.5, 0.1)
  pCO2 = rnorm(nsynth, 280, 10)
  
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
  
  #Diffusion ratio of carbon isotopes 13C/12C
  DIF.ratio = 1.004443
  
  #Isotope ratio constants
  RC.vpdb = 0.011237
  RO.vsmow = 0.0020052
  RO.vpdb = 0.002067181
  
  #Climate Parameters
  PPCQ <- Pa * PfPCQ
  #At least 1mm rain in the carbonate precipitation quarter
  PPCQ = max(PPCQ, 1)
  TmPCQ <- Tma + TmPCQ_min_a
  TmPCQ_K <- TmPCQ + 273.15
  TmOOS <- (Tma * 4 -  TmPCQ) / 3
  
  #Depth to carbonate in meters
  z_m <- z / 100
  
  #Soil temperatures at depth z
  #Assume thermal conductivity = 0.0007 cal / cm2  s  *C, volumetric heat capacity of 0.3 cal / cm2 *C, Quade 2013
  d <- sqrt((2 * 0.0007) / ((2 * 3.1415 / 3.154e7) * 0.3))
  #T is highest (avg summer temps) at t = 0.29, mean at t = 0 (avg spring and fall temps), low at t = 0.79 (avg winter temps)
  t <- ifelse(dryq == 3, 0.29, ifelse(dryq == 1, 0.79, 0))
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
  
  #Respiration rate - Raich 2002, gC/m2day
  R_PCQ_D_m <- 1.25 * exp(0.0545 * TmPCQ) * PPCQ / (127.77 + PPCQ)
  R_PCQ_D_m <- pmax(R_PCQ_D_m * rr, 0.00001)
  R_PCQ_D_m <- R_PCQ_D_m * rr
  #Using mean residual of 50% based on Raich validation data 
  theta = (R_PCQ_D_m * 0.5) ^ 2 / R_PCQ_D_m 
  shp = R_PCQ_D_m / theta
  R_PCQ_D = rgamma(nsynth, shape = shp, scale = theta)
  #Convert to molC/cm^3s
  R_PCQ_D = R_PCQ_D / (12.01 * 100^2)  #molC / cm2 / day
  R_PCQ_S = R_PCQ_D / (24 * 3600)  #molC/ cm2 / s
  R_PCQ_S_0 = R_PCQ_S / (L * pores) # Quade et al. (2007)
  
  #CO2 Diffusion coefficient - based on temp, FAP, tort
  DIFC = FAP * tort * 0.1369 * (TmPCQ_soil_K / 273.15) ^ 1.958
  
  #Water limitation effect of discriminaton, Diefendorf 2010, recalculated by Schubert and Jahren and converted into a correction from water "saturation"
  W_m <- 22.65 - (1.2 * (Pa + 975)) / (27.2 + 0.04 * (Pa + 975))
  W = rnorm(nsynth, W_m, 2.40)
  
  #CO2 effect on discrimination, Schubert and Jahren - bulk above ground tissue
  d13C_pf_pCO2_m <- 28.26 * 0.35 * (pCO2 + 15) / (28.26 + 0.35 * (pCO2 + 15))
  d13C_pf_pCO2 = rnorm(nsynth, d13C_pf_pCO2_m, 0.61)
  
  #Plant discrimination - Assumes bulk above-ground plant tissue approximates d13C values of respiration. Root respiration is depleted compared to SOM and root tissue d13C - Bowling et al. (2007), 
  d13C_R <- d13C_atmCO2 - (d13C_pf_pCO2 - W)
  
  #Convert pCO2 to units mol/cm^3
  pCO2_mcc = pCO2 / (0.08206 * TmPCQ_soil_K * 10^9)
  
  #Soil CO2 C isotopes
  d13C_atmCO2_hat = (d13C_atmCO2 / 1000 + 1) * RC.vpdb / (1 + RC.vpdb * (d13C_atmCO2 / 1000 + 1))
  d13C_R_hat = (d13C_R / 1000 + 1) * RC.vpdb / (1 + RC.vpdb * (d13C_R / 1000 + 1))
  S_z = k ^ 2 * R_PCQ_S_0 / DIFC * (1 - exp(-z / k))
  dC_Soil.num = S_z * DIF.ratio * d13C_R_hat + pCO2_mcc * d13C_atmCO2_hat
  dC_Soil.denom = S_z * (1 - DIF.ratio * d13C_R_hat) + pCO2_mcc * (1 - d13C_atmCO2_hat)
  dC_Soil = (dC_Soil.num / (dC_Soil.denom * RC.vpdb) - 1) * 1000
  R_Soil = (dC_Soil / 1000 + 1) * RC.vpdb
  
  S_z_ppm <- S_z * 0.08206 * TmPCQ_soil_K * 10^9
  
  # Romanek et. al. 1992 Fractionation factor CO2 - calcite
  A_CO2_Carb <- 1 / (1.01198 - 1.2e-4 * TmPCQ_soil) 
  R_Carb <- R_Soil / A_CO2_Carb
  
  #Carb isotope values
  dC_Carb <- (R_Carb / RC.vpdb - 1) * 1000
  
  dat = c(median(dC_Carb), median(S_z_ppm))
  
  return(dat)
}


#Set up respiration ratio optimization
rr = seq(0.01, 1, 0.01)

parms_wq_C.depth = data.frame(rr = rr, rmse = numeric(100), S_z = numeric(100))

## Warm Quarter rr optimization

for(j in 1:nrow(parms_wq_C.depth)){
  opt = data.frame(d13C = numeric(nrow(sites.depth)), S_z = numeric(nrow(sites.depth)))
  for(i in 1: nrow(sites.depth)){
    opt[i,] = sm_optimizer_r.depth(sites.depth$Pa[i], sites.depth$Tma[i], sites.depth$PfWQ[i], sites.depth$TmWQ[i], sites.depth$Ra[i], sites.depth$TmWQ[i], 3, parms_wq_C.depth$rr[j], sites.depth$depth.cm[i])
  }
  
  opt = data.frame(Site = sites.depth$Site, d13C = opt$d13C, S_z = opt$S_z)
  opt = merge.data.frame(opt, sites.depth, by.x = "Site", by.y = "Site", all.x=TRUE)
  
  mse = (opt$d13C - opt$d13C.measured)^2
  mse = mean(mse)
  parms_wq_C.depth$rmse[j] = sqrt(mse)
  parms_wq_C.depth$S_z[j] = mean(opt$S_z)
}

S_z_opt_wq.depth = parms_wq_C.depth[which.min(parms_wq_C.depth$rmse),"S_z"]
rr_opt_wq.depth = parms_wq_C.depth[which.min(parms_wq_C.depth$rmse),"rr"]

S_z_opt_wq.depth
rr_opt_wq.depth

## Dry quarter respiration optimization

parms_dq_C.depth = data.frame(rr = rr, rmse = numeric(100), S_z = numeric(100))

for(j in 1:nrow(parms_dq_C.depth)){
  opt = data.frame(d13C = numeric(nrow(sites.depth)), S_z = numeric(nrow(sites.depth)))
  for(i in 1: nrow(sites.depth)){
    opt[i,] = sm_optimizer_r.depth(sites.depth$Pa[i], sites.depth$Tma[i], sites.depth$PfDQ[i], sites.depth$TmDQ[i], sites.depth$Ra[i], sites.depth$TmWQ[i], sites.depth$DQ[i], parms_wq_C.depth$rr[j], sites.depth$depth.cm[i])
  }
  
  opt = data.frame(Site = sites.depth$Site, d13C = opt$d13C, S_z = opt$S_z)
  opt = merge.data.frame(opt, sites.depth, by.x = "Site", by.y = "Site", all.x=TRUE)
  
  mse = (opt$d13C - opt$d13C.measured)^2
  mse = mean(mse)
  parms_dq_C.depth$rmse[j] = sqrt(mse)
  parms_dq_C.depth$Sz[j] = mean(opt$S_z)
}

S_z_opt_dq.depth = parms_dq_C.depth[which.min(parms_dq_C.depth$rmse),"S_z"]
rr_opt_dq.depth = parms_dq_C.depth[which.min(parms_dq_C.depth$rmse),"rr"]

S_z_opt_dq.depth
rr_opt_dq.depth

#Graph

jpeg("Resp_OPT_depth.jpg", units="in", res=300, width=5, height=4.5)

par(mar=c(4,5,1,1))
plot(parms_wq_C.depth$rmse ~ parms_wq_C.depth$rr, cex = 1.5, cex.axis=1.3, cex.lab=1.3, ylim=c(2,13), type="n", main = "", ylab = expression(paste("RMSE (\u2030)")), xlab = expression(paste("f"[R]," (Fraction of Estimated Respiration)")))
lines(parms_wq_C.depth$rmse ~ parms_wq_C.depth$rr, col="red")
lines(parms_dq_C.depth$rmse ~ parms_dq_C.depth$rr, col="blue")
legend("topright", legend=c("Dry Quarter", "Warm Quarter"), col=c("blue","red"), lty=1)

dev.off()


