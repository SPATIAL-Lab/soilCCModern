## Standalone Optimized models for warm quarter and dry quarter as the PCQ

#number of iterations to calculate error and set a random number generator
nsynth = 10000
set.seed(1)

# Warm quarter model
# Mean annual precipitation (mm / yr)
Pa = 
# Mean annual temperature(deg C)
Tma = 
# Fraction of precipitation in the pedogenic carbonate quarter (as a fraction of Pa)
PfPCQ = 
# Difference between pedogenic carbaonte quarter mean temperature and mean annual temperature
TmPCQ_min_a = 
# Radiation at the top of the atmosphere
Ra = 
# Partial pressure of atmospheric CO2 (ppmV)
pCO2 = 
# Carbon isotope value of atmospheric CO2 (permil)
d13C_atmCO2 = 

carb_model_wq = function(Pa, Tma, PfPCQ, TmPCQ_min_a, Ra, pCO2, d13C_atmCO2){
  
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
  
  #Evaporated soil water optimized
  esw = 0.95
  
  #Seasonal precipitation optimized
  spre = 0.55
  
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
  TmWQ_min_a <- TmPCQ_min_a
  
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
  t <- 0.29
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
  
  # Calculating L (production depth) calculated from aridity index Yang et al. (2016) - approximated from Figure 4a. Assumes production depth is equal to rooting depth - Quade et al. (2007).
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
  
  #Respiration rate
  R_PCQ_D_m <- 1.25 * exp(0.0545 * TmPCQ) * PPCQ / (127.77 + PPCQ)  #Raich 2002, gC/m2day
  R_PCQ_D_m <- R_PCQ_D_m * 0.11
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
  
  # Romanek et. al. 1992 Fractionation factor CO2 - calcite
  A_CO2_Carb <- 1 / (1.01198 - 1.2e-4 * TmPCQ_soil) 
  R_Carb <- R_Soil / A_CO2_Carb
  
  #Precipitation O isotope ratios - WIDB, waterisotopes.org 
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
  dC_Carb <- (R_Carb / RC.vpdb - 1) * 1000
  dO_Carb <- (R_O_Carb / RO.vpdb - 1) * 1000
  
  dat = c(median(dC_Carb), median(dO_Carb), sd(dC_Carb), sd(dO_Carb))
  
  return(dat)
}

wq_isotope_pred = data.frame(d13C = numeric(0), d18O = numeric(0), d13C_sd = numeric(0), d18O_sd = numeric(0))
wq_isotope_pred[1,] = carb_model_wq(Pa, Tma, PfPCQ, TmPCQ_min_a, Ra, pCO2, d13C_atmCO2)
View(wq_isotope_pred)

# Dry quarter model. 
# Mean annual precipitation (mm / yr)
Pa = 
# Mean annual temperature(deg C)
Tma = 
# Fraction of precipitation in the pedogenic carbonate quarter (as a fraction of Pa)
PfPCQ = 
# Difference between pedogenic carbaonte quarter mean temperature and mean annual temperature
TmPCQ_min_a = 
# Radiation at the top of the atmosphere
Ra = 
# Partial pressure of atmospheric CO2 (ppmV)
pCO2 = 
# Carbon isotope value of atmospheric CO2 (permil)
d13C_atmCO2 = 
# Difference between warm quarter mean temperature and mean annual temperature
TmWQ_min_a = 
#DQ dry quarter where is 1 for winter, 2 for spring, 3 for summer and 4 for fall.  
DQ = 

carb_model_dq = function(Pa, Tma, PfPCQ, TmPCQ_min_a, Ra, pCO2, d13C_atmCO2, TmWQ_min_a, DQ){
  
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
  
  #Evaporated soil water optimized
  esw = 0.75
  
  #Seasonal precipitation optimized
  spre = 0.8
  
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
  
  # Calculating L (production depth) calculated from aridity index Yang et al. (2016) - approximated from Figure 4a. Assumes production depth is equal to rooting depth - Quade et al. (2007).
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
  
  #Respiration rate
  R_PCQ_D_m <- 1.25 * exp(0.0545 * TmPCQ) * PPCQ / (127.77 + PPCQ)  #Raich 2002, gC/m2day
  R_PCQ_D_m <- R_PCQ_D_m * 0.19
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
  
  # Romanek et. al. 1992 Fractionation factor CO2 - calcite
  A_CO2_Carb <- 1 / (1.01198 - 1.2e-4 * TmPCQ_soil) 
  R_Carb <- R_Soil / A_CO2_Carb
  
  #Precipitation O isotope ratios - WIDB, waterisotopes.org 
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
  dC_Carb <- (R_Carb / RC.vpdb - 1) * 1000
  dO_Carb <- (R_O_Carb / RO.vpdb - 1) * 1000
  
  dat = c(median(dC_Carb), median(dO_Carb), sd(dC_Carb), sd(dO_Carb))
  
  return(dat)
}

dq_isotope_pred = data.frame(d13C = numeric(0), d18O = numeric(0), d13C_sd = numeric(0), d18O_sd = numeric(0))
dq_isotope_pred[1,] = carb_model_wq(Pa, Tma, PfPCQ, TmPCQ_min_a, Ra, pCO2, d13C_atmCO2, TmWQ_min_a, DQ)
View(dq_isotope_pred)
