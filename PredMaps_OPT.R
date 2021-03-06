library(raster)
library(scales)
library(rgdal)

# Read in rasters - this may make your computer slow
load("rasters.RData")

# Build fuctions for each isotope system and PCQ season
raster_opt_wq_dC = function(Pa, Tma, PfPCQ, TmPCQ_min_a, Ra){
  
  #Pre-industrial pCO2 and carbon isotope value of atmospheric CO2
  d13C_atmCO2 = -6.5
  pCO2 = 280
  
  #Porosity and tortuosity means of 0.35 and 0.7, respectively, and sd of 0.1.
  pores = 0.35
  tort = 0.7
  
  # Solar radiation estimate - Hargreaves and Samani (1982)
  Rs <- Ra * 0.16 * sqrt(12)
  
  #Diffusion ratio of carbon isotopes 13C/12C
  DIF.ratio = 1.004443
  
  #Isotope ratio constants
  RC.vpdb = 0.011237
  RO.vsmow = 0.0020052
  RO.vpdb = 0.002067181
  
  #Climate Parameters
  Pa <- pmax(Pa, 5)
  PPCQ <- Pa * PfPCQ
  #At least 1mm rain in the carbonate precipitation quarter
  PPCQ = pmax(PPCQ, 1)
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
  z = z_min + z_thick / 2 
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
  h <- 0.25 + 0.7 * (PPCQ / 900)
  RH <- h * 100
  
  #Potential Evapotranspiration
  #PET in mm/day, Turc 1961
  #Annual
  PET_A_D <- ifelse (!is.na(Pa) | RH < 50, 0.013 * (Tma / (Tma + 15)) * (23.8856 * Rs + 50)* (1 + ((50 - RH) / 70)), 0.0133 * (Tma / (Tma + 15)) * (23.885 * Rs + 50))
  PET_A_D <- pmax(PET_A_D, 0.01)
  PET_A_A <- PET_A_D * 365
  
  #PCQ
  PET_PCQ_D <- ifelse (!is.na(Pa) | RH < 50, 0.013 * (TmPCQ / (TmPCQ + 15)) * (23.8856 * Rs + 50)* (1 + ((50 - RH) / 70)), 0.0133 * (TmPCQ / (TmPCQ + 15)) * (23.885 * Rs + 50))
  PET_PCQ_D <- pmax(PET_PCQ_D, 0.01)
  PET_PCQ <- PET_PCQ_D * 90
  
  # Calculating L (production depth) calculated from aridity index Yang et al. (2016) - approxiTmaed from Figure 4a. Assumes production depth is equal to rooting depth - Quade et al. (2007).
  AI = PET_A_A / Pa
  L = ifelse(!is.na(Pa) | AI > 1.4, 60, -200 * AI ^ 2 + 250 * AI + 100)
  
  # Production depth relationship with mean production depth and characteristic production depth - Quade et al. (2007)
  k = L / 2 / log(2)
  
  #Actual Evapotranspiration
  #This noise parmeter limits ETA<CMP_mm but allows variation around PET, as observed
  #AET in mm/quarter from Budyko curve
  AET_PCQ = PPCQ * (1 / (sqrt(1 + (1 / ((PET_PCQ / (PPCQ)))) ^ 2))) 
  
  #Free air porosity
  #Scales volumetrically w/ excess precipitation relative to pore space, assume a minimum of 5% volumetric water content
  FAP <- pmin((pores - ((PPCQ - AET_PCQ) / (L * 10 * pores))), pores - 0.05)
  #At least 1% free air porosity
  FAP = pmax(FAP, 0.01)
  
  #Respiration rate
  R_PCQ_D <- 1.25 * exp(0.0545 * TmPCQ) * PPCQ / (127.77 + PPCQ)  #Raich 2002, gC/m2day
  R_PCQ_D <- R_PCQ_D * rr_opt_wq
  #Convert to molC/cm^3s
  R_PCQ_D = R_PCQ_D / (12.01 * 100^2)  #molC / cm2 / day
  R_PCQ_S = R_PCQ_D / (24 * 3600)  #molC/ cm2 / s
  R_PCQ_S_0 = R_PCQ_S / (L * pores) # Quade et al. (2007)
  
  #CO2 Diffusion coefficient - based on temp, FAP, tort
  DIFC = FAP * tort * 0.1369 * (TmPCQ_soil_K / 273.15) ^ 1.958
  
  #Water limitation effect of discriminaton, Diefendorf 2010, recalculated by Schubert and Jahren and converted into a correction from water "saturation"
  W <- 22.65 - (1.2 * (Pa + 975)) / (27.2 + 0.04 * (Pa + 975))
  
  #CO2 effect on discrimination, Schubert and Jahren - bulk above ground tissue
  d13C_pf_pCO2 <- 28.26 * 0.35 * (pCO2 + 15) / (28.26 + 0.35 * (pCO2 + 15))
  
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
  
  #Carb isotope values
  dC_Carb <- (R_Carb / RC.vpdb - 1) * 1000
  
  return(dC_Carb)
}

raster_opt_dq_dC = function(Pa, Tma, PfPCQ, TmPCQ_min_a, Ra, TmWQ_min_a, DQ){
  
  #Pre-industrial pCO2 and carbon isotope value of atmospheric CO2
  d13C_atmCO2 = -6.5
  pCO2 = 280
  
  #Porosity and tortuosity means of 0.35 and 0.7, respectively, and sd of 0.1.
  pores = 0.35
  
  tort = 0.7
  
  # Solar radiation estimate - Hargreaves and Samani (1982)
  Rs <- Ra * 0.16 * sqrt(12)
  
  #Diffusion ratio of carbon isotopes 13C/12C
  DIF.ratio = 1.004443
  
  #Isotope ratio constants
  RC.vpdb = 0.011237
  RO.vsmow = 0.0020052
  RO.vpdb = 0.002067181
  
  #Climate Parameters
  Pa <- pmax(Pa, 5)
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
  z = z_min + z_thick / 2 
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
  h <- 0.25 + 0.7 * (PPCQ / 900)
  RH <- h * 100
  
  #Potential Evapotranspiration
  #PET in mm/day, Turc 1961
  #Annual
  PET_A_D <- ifelse (RH < 50, 0.013 * (Tma / (Tma + 15)) * (23.8856 * Rs + 50)* (1 + ((50 - RH) / 70)), 0.0133 * (Tma / (Tma + 15)) * (23.885 * Rs + 50))
  PET_A_D <- pmax(PET_A_D, 0.01)
  PET_A_A <- PET_A_D * 365
  
  #PCQ
  PET_PCQ_D <- ifelse (RH < 50, 0.013 * (TmPCQ / (TmPCQ + 15)) * (23.8856 * Rs + 50)* (1 + ((50 - RH) / 70)), 0.0133 * (TmPCQ / (TmPCQ + 15)) * (23.885 * Rs + 50))
  PET_PCQ_D <- pmax(PET_PCQ_D, 0.01)
  PET_PCQ <- PET_PCQ_D * 90
  
  # Calculating L (production depth) calculated from aridity index Yang et al. (2016) - approxiTmaed from Figure 4a. Assumes production depth is equal to rooting depth - Quade et al. (2007).
  AI = PET_A_A / Pa
  L = ifelse(AI > 1.4, 60, -200 * AI ^ 2 + 250 * AI + 100)
  
  # Production depth relationship with mean production depth and characteristic production depth - Quade et al. (2007)
  k = L / 2 / log(2)
  
  #Actual Evapotranspiration
  #This noise parmeter limits ETA<CMP_mm but allows variation around PET, as observed
  #AET in mm/quarter from Budyko curve
  AET_PCQ = PPCQ * (1 / (sqrt(1 + (1 / ((PET_PCQ / (PPCQ)))) ^ 2))) 
  
  #Free air porosity
  #Scales volumetrically w/ excess precipitation relative to pore space, assume a minimum of 5% volumetric water content
  FAP <- pmin((pores - ((PPCQ - AET_PCQ) / (L * 10 * pores))), pores - 0.05)
  #At least 1% free air porosity
  FAP = pmax(FAP, 0.01)
  
  #Respiration rate
  R_PCQ_D <- 1.25 * exp(0.0545 * TmPCQ) * PPCQ / (127.77 + PPCQ)  #Raich 2002, gC/m2day
  R_PCQ_D <- R_PCQ_D * rr_opt_dq
  #Convert to molC/cm^3s
  R_PCQ_D = R_PCQ_D / (12.01 * 100^2)  #molC / cm2 / day
  R_PCQ_S = R_PCQ_D / (24 * 3600)  #molC/ cm2 / s
  R_PCQ_S_0 = R_PCQ_S / (L * pores) # Quade et al. (2007)
  
  #CO2 Diffusion coefficient - based on temp, FAP, tort
  DIFC = FAP * tort * 0.1369 * (TmPCQ_soil_K / 273.15) ^ 1.958
  
  #Water limitation effect of discriminaton, Diefendorf 2010, recalculated by Schubert and Jahren and converted into a correction from water "saturation"
  W <- 22.65 - (1.2 * (Pa + 975)) / (27.2 + 0.04 * (Pa + 975))

  #CO2 effect on discrimination, Schubert and Jahren - bulk above ground tissue
  d13C_pf_pCO2 <- 28.26 * 0.35 * (pCO2 + 15) / (28.26 + 0.35 * (pCO2 + 15))
  
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
  
  #Carb isotope values
  dC_Carb <- (R_Carb / RC.vpdb - 1) * 1000
  
  return(dC_Carb)
}

raster_opt_wq_dO = function(Pa, Tma, PfPCQ, TmPCQ_min_a, Ra){
  
  #Porosity and tortuosity means of 0.35 and 0.7, respectively
  pores = 0.35
  tort = 0.7
  
  # Solar radiation estimate - Hargreaves and Samani (1982)
  Rs <- Ra * 0.16 * sqrt(12)
  
  #Isotope ratio constants
  RO.vsmow = 0.0020052
  RO.vpdb = 0.002067181
  
  #Climate Parameters
  Pa <- pmax(Pa, 5)
  PPCQ <- Pa * PfPCQ
  #At least 1mm rain in the carbonate precipitation quarter
  PPCQ = pmax(PPCQ, 1)
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
  z = z_min + z_thick / 2 
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
  h <- 0.25 + 0.7 * (PPCQ / 900)
  RH <- h * 100
  
  #Potential Evapotranspiration
  #PET in mm/day, Turc 1961
  #Annual
  PET_A_D <- ifelse (RH < 50, 0.013 * (Tma / (Tma + 15)) * (23.8856 * Rs + 50)* (1 + ((50 - RH) / 70)), 0.0133 * (Tma / (Tma + 15)) * (23.885 * Rs + 50))
  PET_A_D <- pmax(PET_A_D, 0.01)
  PET_A_A <- PET_A_D * 365
  
  #PCQ
  PET_PCQ_D <- ifelse (RH < 50, 0.013 * (TmPCQ / (TmPCQ + 15)) * (23.8856 * Rs + 50)* (1 + ((50 - RH) / 70)), 0.0133 * (TmPCQ / (TmPCQ + 15)) * (23.885 * Rs + 50))
  PET_PCQ_D <- pmax(PET_PCQ_D, 0.01)
  PET_PCQ <- PET_PCQ_D * 90
  
  # Calculating L (production depth) calculated from aridity index Yang et al. (2016) - approxiTmaed from Figure 4a. Assumes production depth is equal to rooting depth - Quade et al. (2007).
  AI = PET_A_A / Pa
  L = ifelse(AI > 1.4, 60, -200 * AI ^ 2 + 250 * AI + 100)
  
  # Production depth relationship with mean production depth and characteristic production depth - Quade et al. (2007)
  k = L / 2 / log(2)
  
  #Actual Evapotranspiration
  #AET in mm/quarter from Budyko curve
  AET_PCQ = PPCQ * (1 / (sqrt(1 + (1 / ((PET_PCQ / (PPCQ)))) ^ 2))) 
  
  #Free air porosity
  #Scales volumetrically w/ excess precipitation relative to pore space, assume a minimum of 5% volumetric water content
  FAP <- pmin((pores - ((PPCQ - AET_PCQ) / (L * 10 * pores))), pores - 0.05)
  #At least 1% free air porosity
  FAP = pmax(FAP, 0.01)
  
  #Precipitation O isotope ratios - WIDB, waterisotopes.org 
  dO_P_OOS <- -15.57 + 0.60 * TmOOS
  dO_P_PCQ <- -15.57 + 0.60 * TmPCQ
  dO_P_soil <- (dO_P_PCQ * PfPCQ + dO_P_OOS * (1 - spre_opt_wq) * (1 - PfPCQ)) / (PfPCQ + (1 - spre_opt_wq) * (1 - PfPCQ))
  R_O_P_soil = (dO_P_soil / 1000 + 1) * RO.vsmow
  R_O_P_PCQ = (dO_P_PCQ / 1000 + 1) * RO.vsmow
  
  #Atmospheric vapor O isotope ratios
  A_atmH2O_P <- 2.71828 ^ ((5.9702e6 / TmPCQ_K ^ 2 - 3.2801e4 / TmPCQ_K + 52.227) / 1000)
  R_O_atmH2O <- R_O_P_PCQ / A_atmH2O_P
  dO_atmH2O <- (R_O_atmH2O / RO.vsmow - 1) * 1000
  
  #Soil water evaporation
  #evap is 6%+/-4% of total ET globally Good et. al. 2010
  E = 0.06 * AET_PCQ
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
  R_O_soil = R_O_soil * esw_opt_wq + R_O_P_soil * (1 - esw_opt_wq)  
  dO_soil <- (R_O_soil/RO.vsmow - 1) * 1000
  
  #Soil carbonate O isotope fractionation - Kim and O'neal 1997
  A_H20_Carb <- 2.71828 ^ ((1.803e4 / TmPCQ_soil_K - 32.42) / 1000)
  R_O_Carb <- R_O_soil * A_H20_Carb
  
  #Carb isotope values
  dO_Carb <- (R_O_Carb / RO.vpdb - 1) * 1000
  
  return(dO_Carb)
}

raster_opt_dq_dO = function(Pa, Tma, PfPCQ, TmPCQ_min_a, Ra, TmWQ_min_a, DQ){
  
  #Porosity and tortuosity means of 0.35 and 0.7, respectively
  pores = 0.35
  tort = 0.7
  
  # Solar radiation estimate - Hargreaves and Samani (1982)
  Rs <- Ra * 0.16 * sqrt(12)
  
  #Isotope ratio constants
  RO.vsmow = 0.0020052
  RO.vpdb = 0.002067181
  
  #Climate Parameters
  Pa <- pmax(Pa, 5)
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
  z = z_min + z_thick / 2 
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
  h <- 0.25 + 0.7 * (PPCQ / 900)
  RH <- h * 100
  
  #Potential Evapotranspiration
  #PET in mm/day, Turc 1961
  #Annual
  PET_A_D <- ifelse (RH < 50, 0.013 * (Tma / (Tma + 15)) * (23.8856 * Rs + 50)* (1 + ((50 - RH) / 70)), 0.0133 * (Tma / (Tma + 15)) * (23.885 * Rs + 50))
  PET_A_D <- pmax(PET_A_D, 0.01)
  PET_A_A <- PET_A_D * 365
  
  #PCQ
  PET_PCQ_D <- ifelse (RH < 50, 0.013 * (TmPCQ / (TmPCQ + 15)) * (23.8856 * Rs + 50)* (1 + ((50 - RH) / 70)), 0.0133 * (TmPCQ / (TmPCQ + 15)) * (23.885 * Rs + 50))
  PET_PCQ_D <- pmax(PET_PCQ_D, 0.01)
  PET_PCQ <- PET_PCQ_D * 90
  
  # Calculating L (production depth) calculated from aridity index Yang et al. (2016) - approxiTmaed from Figure 4a. Assumes production depth is equal to rooting depth - Quade et al. (2007).
  AI = PET_A_A / Pa
  L = ifelse(AI > 1.4, 60, -200 * AI ^ 2 + 250 * AI + 100)
  
  # Production depth relationship with mean production depth and characteristic production depth - Quade et al. (2007)
  k = L / 2 / log(2)
  
  #Actual Evapotranspiration
  #AET in mm/quarter from Budyko curve
  AET_PCQ = PPCQ * (1 / (sqrt(1 + (1 / ((PET_PCQ / (PPCQ)))) ^ 2))) 
  
  #Free air porosity
  #Scales volumetrically w/ excess precipitation relative to pore space, assume a minimum of 5% volumetric water content
  FAP <- pmin((pores - ((PPCQ - AET_PCQ) / (L * 10 * pores))), pores - 0.05)
  #At least 1% free air porosity
  FAP = pmax(FAP, 0.01)
  
  #Precipitation O isotope ratios - WIDB, waterisotopes.org 
  dO_P_OOS <- -15.57 + 0.60 * TmOOS
  dO_P_PCQ <- -15.57 + 0.60 * TmPCQ
  dO_P_soil <- (dO_P_PCQ * PfPCQ + dO_P_OOS * (1 - spre_opt_dq) * (1 - PfPCQ)) / (PfPCQ + (1 - spre_opt_dq) * (1 - PfPCQ))
  R_O_P_soil = (dO_P_soil / 1000 + 1) * RO.vsmow
  R_O_P_PCQ = (dO_P_PCQ / 1000 + 1) * RO.vsmow
  
  #Atmospheric vapor O isotope ratios
  A_atmH2O_P <- 2.71828 ^ ((5.9702e6 / TmPCQ_K ^ 2 - 3.2801e4 / TmPCQ_K + 52.227) / 1000)
  R_O_atmH2O <- R_O_P_PCQ / A_atmH2O_P
  dO_atmH2O <- (R_O_atmH2O / RO.vsmow - 1) * 1000
  
  #Soil water evaporation
  #evap is 6%+/-4% of total ET globally Good et. al. 2010
  E = 0.06 * AET_PCQ
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
  R_O_soil = R_O_soil * esw_opt_dq + R_O_P_soil * (1 - esw_opt_dq)  
  dO_soil <- (R_O_soil/RO.vsmow - 1) * 1000
  
  #Soil carbonate O isotope fractionation - Kim and O'neal 1997
  A_H20_Carb <- 2.71828 ^ ((1.803e4 / TmPCQ_soil_K - 32.42) / 1000)
  R_O_Carb <- R_O_soil * A_H20_Carb
  
  #Carb isotope values
  dO_Carb <- (R_O_Carb / RO.vpdb - 1) * 1000
  
  return(dO_Carb)
}

#Vectorize the functions so we can use overlay
raster_opt_wq_dC <- Vectorize(raster_opt_wq_dC)
raster_opt_dq_dC <- Vectorize(raster_opt_dq_dC)
raster_opt_wq_dO <- Vectorize(raster_opt_wq_dO)
raster_opt_dq_dO <- Vectorize(raster_opt_dq_dO)


#Use overlay to do complex raster calcs
dC_Carb_map_wq <- overlay(Pa, Tma, PfWQ, TmWQ_min_a, Ra, fun=raster_opt_wq_dC)
dC_Carb_map_dq <- overlay(Pa, Tma, PfDQ, TmDQ_min_a, Ra, TmWQ_min_a, DQ, fun=raster_opt_dq_dC)
dO_Carb_map_wq <- overlay(Pa, Tma, PfWQ, TmWQ_min_a, Ra, fun=raster_opt_wq_dO)
dO_Carb_map_dq <- overlay(Pa, Tma, PfDQ, TmDQ_min_a, Ra, TmWQ_min_a, DQ, fun=raster_opt_dq_dO)

dO_Carb_map_diff <- dO_Carb_map_wq - dO_Carb_map_dq
dC_Carb_map_diff <- dC_Carb_map_wq - dC_Carb_map_dq

# Raster of plant discrimination / difference between d13Cplant and d13Ccarb
pCO2 <- 280
d13C_atmCO2 <- -6.5

#Water limitation of discriminaton raster, Diefendorf
W <- 22.65 - (1.2 * (Pa + 975)) / (27.2 + 0.04 * (Pa + 975))

#CO2 effect on discrimination raster, Schubert and Jahren
d13C_pf_pCO2 <- 28.26 * 0.35 * (pCO2 + 15) / (28.26 + 0.35 * (pCO2 + 15))

#Plant discrimination raster
d13C_R <- d13C_atmCO2 - (d13C_pf_pCO2 - W)
plant_carb_diff_map <- dC_Carb_map_wq - d13C_R

#Raster of difference between mean precipitation and carbonate o isotope ratio
dO_P_map <- -15.57 + 0.60 * Tma
RO.vpdb <- 0.002067181
RO.vsmow <- 0.0020052
RO_P_map <- (dO_P_map / 1000 + 1) * RO.vsmow
A_H20_Carb <- 2.71828 ^ ((1.803e4 / (Tma + 273.15) - 32.42) / 1000)
R_O_Carb <- RO_P_map * A_H20_Carb
dO_Carb_a <- (R_O_Carb / RO.vpdb - 1) * 1000

dO_diff_map <- dO_Carb_map_wq - dO_Carb_a

# Delete C4 < 10 % veg cover
C4_10 <- C4
C4_10[C4 < 10] <- NA
plot(C4_10)

#jpeg("PredMaps_OPT.jpg", units="in", res=300, width=9.8, height=7)

layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow=TRUE), heights=c(1,1,1,1), widths=c(1,1,1,1))

map_cols_C4 <- colorRampPalette(c(alpha("white", 0), alpha("grey", 0.3), alpha("grey", 0.4), alpha("grey", 0.5)))

par(mar=c(1,0,1,4))
plot(dO_Carb_map_wq, xaxt = 'n', yaxt = 'n')
title("")
mtext(expression(paste(delta^{18}, "O"[wq] ," (\u2030)")), 4, line=1.1, cex=1)
text(-170, 80 ,"a", cex=1.5)
plot(C4_10, col=grey(c(0.4,0.3,0.2), alpha=0.4), add = TRUE, legend=F)

par(mar=c(1,0,1,4))
plot(dC_Carb_map_wq, xaxt = 'n', yaxt = 'n')
title("")
mtext(expression(paste(delta^{13}, "C"[wq] ," (\u2030)")), 4, line=1.1, cex=1)
text(-170, 80 ,"b", cex=1.5)
plot(C4_10, col=grey(c(0.4,0.3,0.2), alpha=0.4), add = TRUE, legend=F)

par(mar=c(1,0,1,4))
plot(dO_Carb_map_diff, xaxt = 'n', yaxt = 'n')
title("")
mtext(expression(paste(Delta^{18}, "O"[wq - dq] ,"(\u2030)")), 4, line=1.1, cex=1)
text(-170, 80 ,"c", cex=1.5)
plot(C4_10, col=grey(c(0.4,0.3,0.2), alpha=0.4), add = TRUE, legend=F)

par(mar=c(1,0,1,4))
plot(dC_Carb_map_diff, xaxt = 'n', yaxt = 'n')
title("")
mtext(expression(paste(Delta^{13}, "C"[wq - dq] ," (\u2030)")), 4, line=1.1, cex=1)
text(-170, 80 ,"d", cex=1.5)
plot(C4_10, col=grey(c(0.4,0.3,0.2), alpha=0.4), add = TRUE, legend=F)

par(mar=c(1,0,1,4))
plot(dO_diff_map, xaxt = 'n', yaxt = 'n')
title("")
mtext(expression(paste(Delta^{18}, "O"[wq - Pa] ," (\u2030)")), 4, line=1.1, cex=1)
text(-170, 80 ,"e", cex=1.5)
plot(C4_10, col=grey(c(0.4,0.3,0.2), alpha=0.4), add = TRUE, legend=F)

par(mar=c(1,0,1,4))
plot(plant_carb_diff_map, xaxt = 'n', yaxt = 'n')
title("")
mtext(expression(paste(Delta^{13}, "C"[wq - plant] ," (\u2030)")), 4, line=1.1, cex=1)
text(-170, 80 ,"f", cex=1.5)
plot(C4_10, col=grey(c(0.4,0.3,0.2), alpha=0.4), add = TRUE, legend=F)

#dev.off()
