setwd("C:/Users/femal/Dropbox/worldclim")

#Prepare gridded climate data
library(raster)
library(rgdal)

#read in worldclim 2.5min mean temperature grids, downloaded 9-18-18

files = list.files(pattern="wc2.0_2.5m_tavg*")
mtemp = stack(files)
mat = mean(mtemp)
writeRaster(mat, "wc2.0_2.5m_tavg_ma.tif")
plot(mat)

#now make quarterly seasonal averages
qtemp = stack(mean(subset(mtemp, c(12,1,2))), 
              mean(subset(mtemp, 3:5)), 
              mean(subset(mtemp, 6:8)), 
              mean(subset(mtemp, 9:11)))

#find the hottest quarter
hotq = which.max(qtemp)
plot(hotq)

#read in worldclim 2.5min precip grids, downloaded 9-18-18
files = list.files(pattern="wc2.0_2.5m_prec*")
mprec = stack(files)
map = sum(mprec)
writeRaster(map, "wc2.0_2.5m_prec_ma.tif", overwrite=TRUE)
plot(map)

#now make quarterly seasonal sums
qprec = stack(sum(subset(mprec, c(12,1,2))), 
              sum(subset(mprec, 3:5)), 
              sum(subset(mprec, 6:8)), 
              sum(subset(mprec, 9:11)))

#find the driest quarter
dryq = which.min(qprec)
plot(dryq)

#now pull hot quarter temps into a single raster
hqt = max(qtemp)
plot(hqt)

#this is what that would look like as tseas
hqt.offset = hqt - mat
plot(hqt.offset)

#and get hot quarter precip as well
#first pull the values from the raster into a matrix
qprec.mat = matrix(c(qprec[[1]][], qprec[[2]][], qprec[[3]][], qprec[[4]][]), length(qprec[[1]]), 4)
#now set up space in a raster to store results
hqp = hotq
#now use the magic of cbind to select the matrix colums for the hot quarter and combine
hqp[] = qprec.mat[cbind(1:nrow(qprec.mat), hotq[])]
plot(hqp)

#then this is what we want as pseas for the model
hqp.frac = hqp/map
plot(hqp.frac)

#now get dry quarter temp, using the same method
qtemp.mat = matrix(c(qtemp[[1]][], qtemp[[2]][], qtemp[[3]][], qtemp[[4]][]), length(qtemp[[1]]), 4)
dqt = dryq
dqt[] = qtemp.mat[cbind(1:nrow(qtemp.mat), dryq[])]
plot(dqt)

dqt.offset = dqt - mat
plot(dqt.offset)

#and lastly the dry quarter precip
dqp = min(qprec)
plot(dqp)

dqp.frac = dqp/map
plot(dqp.frac)

### Get rid of some of the rasters we don't need
qprec <- NULL
qtemp <- NULL
dqt <- NULL
qtemp.mat <- NULL
qprec.mat <- NULL
hqp <- NULL
hqt<- NULL
dryq <- NULL
hotq <- NULL
dqp <- NULL

#####That's all the grid processing - next bit needs only be run once after data are updated

#read in the compiled data
setwd("~/GitHub/soilCCModern/")
library(xlsx)

#extract relevant values at sites
sites = read.xlsx("modern_comparison.xlsx", sheetIndex = 1)
coords = matrix(c(sites$Lon, sites$Lat),nrow(sites),2)
sites$mat.wc = extract(mat, coords)
sites$map.wc = extract(map, coords)
sites$hqt.offset = extract(hqt.offset, coords)
sites$hqp.frac = extract(hqp.frac, coords)
sites$dqt.offset = extract(dqt.offset, coords)
sites$dqp.frac = extract(dqp.frac, coords)

write.csv(sites, "valsites.csv")

#This does the same for a screened subset of sites
sites = read.xlsx("modern_comparison.xlsx", sheetIndex = 2)
coords = matrix(c(sites$Lon, sites$Lat),nrow(sites),2)
sites$mat.wc = extract(mat, coords)
sites$map.wc = extract(map, coords)
sites$hqt.offset = extract(hqt.offset, coords)
sites$hqp.frac = extract(hqp.frac, coords)
sites$dqt.offset = extract(dqt.offset, coords)
sites$dqp.frac = extract(dqp.frac, coords)

write.csv(sites, "valsites_sel.csv")

#some code to plot up values for some sites
plot(sites$mat.wc, sites$map.wc)
plot(sites$MAT, sites$mat.wc)
abline(0,1)

plot(sites$MAP, sites$map.wc)
abline(0,1)

stack(mat)
## Create prediction maps for oxygen and carbon isotopes of pedogenic carbonates to test against
# Create raster for estimating radiation values at each latitude

Ra <- raster(ncol=8640, nrow=4320, res=c(0.04166667, 0.04166667), xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

for (i in 1:ncol(Ra)){
  
  lat = 0.0416667 * (i - 1 - 4320/2)
  Ra[i,] = 34.5 * cos(abs(lat * 0.01745))
  
}

plot(Ra)


raster_opt_hq_dC = function(map, mat, hqp.frac, hqt.offset, Ra) {
  
  MAP = map
  MAT = mat
  T_seas = hqt.offset
  P_seas = hqp.frac
  
  pCO2 = 280
  deltaA = -6.5
  pores = 0.4
  tort = 0.7
  
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
  
  #Soil temps - T is highest at t = 0.3, mean at t = 0.05, low at t = 0.8
  
  t <- 0.3
  d <- sqrt((2*0.0007)/((2*3.1415/3.154e7)*0.3))
  T_soil <- MAT + (T_seas * sin(2 * 3.1415 * t - z/d) / exp(z/d)) 
  T_soil_K <- T_soil + 273.15
  
  #Convert pCO2 to units mol/cm^3
  pCO2_mcc = pCO2 / (0.08206 * T_soil_K * 10^9)  #mol/cm^3
  
  #Relative humidity, now beta distribution
  h <- 0.25 + 0.7 * (CQP / 900)
  RH <- h * 100
  
  #Depth to carbonate, now gamma dist
  z_mean <- MAP * 0.0925 + 13.4  #top of Bk equation based on Retallack 2005 data
  z_thick = abs(CMP_mm - MAP/12) * 0.74 + 17.3 #thickness of Bk, CMP as proxy for seasonality 
  z = z_mean + z_thick/2 #find middle of Bk in cm
  z <- pmin(z, 100) 
  
  #Respiration rate, now gamma dist
  R_day_m <- 1.25 * exp(0.0545 * CQT) * CMP_cm / (4.259 + CMP_cm)  #Raich 2002, gC/m2day
  R_day <- R_day_m * rr_opt_hq # Hot quarter respiration rate optimization - 39% of estimated
  R_day = R_day / (12.01 * 100^2)  #molC/cm2day
  R_sec <- R_day / (24 * 3600)  #molC/cm2s
  R_sec = R_sec / L / pores #molC/cm3s
  
  #
  
  #Potential ET
  ETP_D <- ifelse (RH < 50, 0.013 * (CQT / (CQT + 15)) * (23.8856 * Rs + 50) * (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50))
  ETP_M <- ETP_D * 30  #mm/month
  
  #Actual ET
  ETA_var = 1 #This noise parmeter limits ETA<CMP_mm but allows variation around ETP, as observed
  ETA = CMP_mm*3 * (1 / (sqrt(1 + (1 / ((ETP_M / (CMP_mm*3)) * ETA_var)) ^ 2))) #AET in mm/month from Budyko curve
  #here scaled eta to quarter precip, assuming potential carry-over
  
  #Free air porosity
  #Have updated, now scales volumetrically w/ excess precipitation relative to pore space
  FAP <- pmin((pores - (CMP_mm - ETA)/(L*10*pores)), pores - 0.05)
  FAP = pmax(FAP,0.01) #dimensionless
  
  #CO2 Diffusion coefficients
  DIFC = FAP * tort * 0.1369 * (T_soil_K / 273.15) ^ 1.958
  
  #Water limitation of discriminaton, Diefendorf
  W <- 22.65 - (1.2 * (MAP + 975)) / (27.2 + 0.04 * (MAP + 975))
  
  #CO2 effect on discrimination, Schubert
  deltaP_pCO2 <- 28.26 * 0.35 * (pCO2 + 15) / (28.26 + 0.35 * (pCO2 + 15))
  
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
  
  return(dC_Carb)
}

raster_opt_dq_dC = function(map, mat, dqp.frac, dqt.offset, hqt.offset, Ra) {
  
  MAP = map
  MAT = mat
  T_seas_hq = hqt.offset
  T_seas = dqt.offset
  P_seas = dqp.frac
  
  pCO2 = 280
  deltaA = -6.5
  pores = 0.4
  tort = 0.7
  
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
  
  #Soil temps - T is highest at t = 0.3, mean at t = 0.05, low at t = 0.8
  
  t <- ifelse(is.na(T_seas) | T_seas_hq == T_seas, 0.3, ifelse(CQT < MAT - 3, 0.8, 0.05))
  d <- sqrt((2*0.0007)/((2*3.1415/3.154e7)*0.3))
  T_soil <- MAT + (T_seas_hq * sin(2 * 3.1415 * t - z/d) / exp(z/d)) 
  T_soil_K <- T_soil + 273.15
  
  #Convert pCO2 to units mol/cm^3
  pCO2_mcc = pCO2 / (0.08206 * T_soil_K * 10^9)  #mol/cm^3
  
  #Relative humidity, now beta distribution
  h <- 0.25 + 0.7 * (CQP / 900)
  RH <- h * 100
  
  #Depth to carbonate, now gamma dist
  z_mean <- MAP*0.0925 + 13.4  #top of Bk equation based on Retallack 2005 data
  z_thick = abs(CMP_mm - MAP/12) * 0.74 + 17.3 #thickness of Bk, CMP as proxy for seasonality 
  z = z_mean + z_thick/2 #find middle of Bk in cm
  z <- pmin(z, 100) 
  
  #Respiration rate, now gamma dist
  R_day_m <- 1.25 * exp(0.0545 * CQT) * CMP_cm / (4.259 + CMP_cm)  #Raich 2002, gC/m2day
  R_day <- R_day_m * rr_opt_hq # Hot quarter respiration rate optimization - 39% of estimated
  R_day = R_day / (12.01 * 100^2)  #molC/cm2day
  R_sec <- R_day / (24 * 3600)  #molC/cm2s
  R_sec = R_sec / L / pores #molC/cm3s
  
  #Potential ET
  ETP_D <- ifelse (RH < 50, 0.013 * (CQT / (CQT + 15)) * (23.8856 * Rs + 50) * (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50))
  ETP_M <- ETP_D * 30  #mm/month
  
  #Actual ET
  ETA_var = 1 #This noise parmeter limits ETA<CMP_mm but allows variation around ETP, as observed
  ETA = CMP_mm*3 * (1 / (sqrt(1 + (1 / ((ETP_M / (CMP_mm*3)) * ETA_var)) ^ 2))) #AET in mm/month from Budyko curve
  #here scaled eta to quarter precip, assuming potential carry-over
  
  #Free air porosity
  #Have updated, now scales volumetrically w/ excess precipitation relative to pore space
  FAP <- pmin((pores - (CMP_mm - ETA)/(L*10*pores)), pores - 0.05)
  FAP = pmax(FAP,0.01) #dimensionless
  
  #CO2 Diffusion coefficients
  DIFC = FAP * tort * 0.1369 * (T_soil_K / 273.15) ^ 1.958
  
  #Water limitation of discriminaton, Diefendorf
  W <- 22.65 - (1.2 * (MAP + 975)) / (27.2 + 0.04 * (MAP + 975))
  
  #CO2 effect on discrimination, Schubert
  deltaP_pCO2 <- 28.26 * 0.35 * (pCO2 + 15) / (28.26 + 0.35 * (pCO2 + 15))
  
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
  
  return(dC_Carb)
}

raster_opt_hq_dO = function(map, mat, hqp.frac, hqt.offset){
  
  MAT = mat
  T_seas = hqt.offset
  CQT <- MAT + T_seas
  CQT_K <- CQT + 273
  MAP = map
  P_seas = hqp.frac
  CQP <- MAP * P_seas
  CQP = max(CQP, 3)
  CMP_mm <- CQP / 3
  CMP_cm <- CQP / 30
  RO.vsmow = 0.0020052
  RO.vpdb = 0.002067181
  
  #Soil temps - T is highest at t = 0.3, mean at t = 0.05, low at t = 0.8
  z_mean <- MAP*0.0925 + 13.4  #top of Bk equation based on Retallack 2005 data
  z_thick = abs(CMP_mm - MAP/12) * 0.74 + 17.3 #thickness of Bk, CMP as proxy for seasonality 
  z = z_mean + z_thick/2 #find middle of Bk in cm
  z <- pmin(z, 100)
  
  t <- 0.3
  d <- sqrt((2*0.0007)/((2*3.1415/3.154e7)*0.3))
  T_soil <- MAT + (T_seas * sin(2 * 3.1415 * t - z/d) / exp(z/d)) 
  T_soil_K <- T_soil + 273.15
  # Precipitation O isotope ratios 
  dO_P <- -14.76 + 0.56 * MAT
  R_O_P = (dO_P / 1000 + 1) * RO.vsmow
  
  #Soil carbonate O isotopes
  A_O <- 2.71828 ^ ((2.78e6 / T_soil_K ^ 2 - 2.89) / 1000)
  R_O_Carb <- R_O_P * A_O
  dO_Carb <- (R_O_Carb / RO.vpdb - 1) * 1000
  
  return(dO_Carb)
  
}

raster_opt_dq_dO = function(map, mat, dqp.frac, dqt.offset, hqt.offset){
  
  MAT = mat
  T_seas = dqt.offset
  T_seas_hq = hqt.offset
  MAP = map
  P_seas = dqp.frac
  CQT <- MAT + T_seas
  CQT_K <- CQT + 273
  CQP <- MAP * P_seas
  CQP = max(CQP, 3)
  CMP_mm <- CQP / 3
  CMP_cm <- CQP / 30
  RO.vsmow = 0.0020052
  RO.vpdb = 0.002067181
  
  #Soil temps - T is highest at t = 0.3, mean at t = 0.05, low at t = 0.8
  z_mean <- MAP * 0.0925 + 13.4  #top of Bk equation based on Retallack 2005 data
  z_thick = abs(CMP_mm - MAP/12) * 0.74 + 17.3 #thickness of Bk, CMP as proxy for seasonality 
  z = z_mean + z_thick/2 #find middle of Bk in cm
  z <- pmin(z, 100)
  
  t <- ifelse(is.na(T_seas) | T_seas_hq == T_seas, 0.3, ifelse(CQT < MAT - 3, 0.8, 0.05))

  d <- sqrt((2*0.0007)/((2*3.1415/3.154e7)*0.3))
  T_soil <- MAT + (T_seas_hq * sin(2 * 3.1415 * t - z/d) / exp(z/d)) 
  T_soil_K <- T_soil + 273.15
  # Precipitation O isotope ratios 
  dO_P <- -14.76 + 0.56 * MAT
  R_O_P = (dO_P / 1000 + 1) * RO.vsmow
  
  #Soil carbonate O isotopes
  A_O <- 2.71828 ^ ((2.78e6 / T_soil_K ^ 2 - 2.89) / 1000)
  R_O_Carb <- R_O_P * A_O
  dO_Carb <- (R_O_Carb / RO.vpdb - 1) * 1000
  
  
  
  return(dO_Carb)
  
}

raster_opt_hq_dO <- Vectorize(raster_opt_hq_dO)
raster_opt_dq_dO <- Vectorize(raster_opt_dq_dO)
raster_opt_hq_dC <- Vectorize(raster_opt_hq_dC)
raster_opt_dq_dC <- Vectorize(raster_opt_dq_dC)

dO_Carb_map_hq <- overlay(map, mat, hqp.frac, hqt.offset, fun=raster_opt_hq_dO)
dC_Carb_map_hq <- overlay(map, mat, hqp.frac, hqt.offset, Ra, fun=raster_opt_hq_dC)
dO_Carb_map_dq <- overlay(map, mat, dqp.frac, dqt.offset, hqt.offset, fun=raster_opt_dq_dO)
dC_Carb_map_dq <- overlay(map, mat, dqp.frac, dqt.offset, hqt.offset, Ra, fun=raster_opt_dq_dC)

jpeg("OPT_MAPS.jpg", units="in", res=300, width=8.5, height=11)

plot(dO_Carb_map_hq)
plot(dC_Carb_map_hq)
plot(dO_Carb_map_dq)
plot(dC_Carb_map_dq)

dev.off()

save.image(file=".RData")

