###Code for conducting validation of soil carbonate model

#Prepare gridded climate data

library(raster)
library(rgdal)
library(readxl)

#read in worldclim 2.5min mean temperature grids, downloaded 9-18-18

setwd("C:/Users/femal/Dropbox/worldclim/")

files = list.files(pattern="wc2.0_2.5m_tavg*")
Tm = stack(files)
Tma = mean(mtemp)
writeRaster(Tma, "wc2.0_2.5m_Tma.tif")

#now make quarterly seasonal averages
TmQ = stack(mean(subset(Tm, c(12,1,2))), 
              mean(subset(Tm, 3:5)), 
              mean(subset(Tm, 6:8)), 
              mean(subset(Tm, 9:11)))

TmQ_djf <- TmQ[[1]]
TmQ_mam <- TmQ[[2]]
TmQ_jja <- TmQ[[3]]
TmQ_son <- TmQ[[4]]

writeRaster(TmQ_djf, "TmQ_djf.tif")
writeRaster(TmQ_mam, "TmQ_mam.tif")
writeRaster(TmQ_jja, "TmQ_jja.tif")
writeRaster(TmQ_son, "TmQ_son.tif")

#find the warmest quarter
WQ = which.max(TmQ)

#read in worldclim 2.5min precip grids, downloaded 9-18-18
files = list.files(pattern="wc2.0_2.5m_prec*")
P = stack(files)
Pa = sum(P)
writeRaster(Pa, "wc2.0_2.5m_Pa.tif", overwrite=TRUE)

#now make quarterly seasonal sums
PQ = stack(sum(subset(P, c(12,1,2))), 
              sum(subset(P, 3:5)), 
              sum(subset(P, 6:8)), 
              sum(subset(P, 9:11)))

#find the driest quarter
DQ = which.min(PQ)

#now pull warm quarter temps into a single raster
TmWQ = max(Tmq)

#this is what that would look like as Tmq_min_a
TmWQ_min_a = TmWQ - Tma

#and get warm quarter precip as well
#first pull the values from the raster into a matrix
PQ.mat = matrix(c(PQ[[1]][], PQ[[2]][], PQ[[3]][], PQ[[4]][]), length(PQ[[1]]), 4)
#now set up space in a raster to store results
PWQ = WQ
#now use the magic of cbind to select the matrix colums for the hot quarter and combine
PWQ[] = PQ.mat[cbind(1:nrow(Pq.mat), WQ[])]

#then this is what we want as pseas for the model
PfWQ = PWQ/Pa

#now get dry quarter temp, using the same method
TmQ.mat = matrix(c(TmQ[[1]][], TmQ[[2]][], TmQ[[3]][], TmQ[[4]][]), length(TmQ[[1]]), 4)
TmDQ = DQ
TmDQ[] = TmQ.mat[cbind(1:nrow(TmQ.mat), DQ[])]

TmDQ_min_a = TmDQ - Tma

#and lastly the dry quarter precip
PDQ = min(PQ)

PfDQ = PDQ/Pa

### To remove raster objects once done with them (they are large and will slow down other processes)
rm(list=ls()[grep("Tm",ls())])
rm(list=ls()[grep("P",ls())])

### Record which quarter is dry quarter. Change southern hemisphere sites to jja = winter, djf = summer etc...

DQ_S_ex <- extent(c(-180, 180, -90, 0))
DQ_S <- crop(DQ, DQ_S_ex)
DQ_N_ex <- extent(c(-180, 180, 0, 90))
DQ_N <- crop(DQ, DQ_N_ex)

DQ_S[DQ == 1] <- 5
DQ_S[DQ == 3] <- 1
DQ_S[DQ_S == 5] <- 3
DQ_S[DQ == 2] <- 5
DQ_S[DQ == 4] <- 2
DQ_S[DQ_S == 5] <- 4

DQ <- merge(DQ_N, DQ_S)

# 1 is winter, 2 is spring, 3 is summer, 4 is fall

# Create raster for estimating avg radiation values at each latitude based on 34.5 Mj / m2 / d at equator and cos(latitude) - rough estimate

Ra <- raster(ncol=8640, nrow=4320, res=c(0.04166667, 0.04166667), xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

for (i in 1:ncol(Ra)){
  
  lat = 0.0416667 * (i - 1 - 4320/2)
  Ra[i,] = 34.5 * cos(abs(lat * 0.01745))
  
}

## Read raster of C4 (as a percentage of veg cover) distribution worldwide - Still et al. (2003)

C4 <- raster("c4_percent_1d.asc")
C4 <- max(C4$c4_percent_1d, 0)
writeRaster(C4, "C4.tif")

## If you have already done these raster calcs, read them in
setwd("C:/Users/femal/Dropbox/worldclim/")
PfWQ <- raster("wc2.0_2.5m_PfWQ.tif")
PfDQ <- raster("wc2.0_2.5m_PfDQ.tif")
TmWQ_min_a <- raster("wc2.0_2.5m_TmWQ_min_a.tif")
TmDQ_min_a <- raster("wc2.0_2.5m_TmDQ_min_a.tif")
Pa <- raster("wc2.0_2.5m_Pa.tif")
Tma <- raster("wc2.0_2.5m_Tma.tif")
DQ <- raster("DQ.tif")
Ra <- raster("Ra.tif")
C4 <- raster("C4.tif")

#####That's all the grid processing - next bit needs only be run once after data are updated

#read in the compiled data
setwd("C:/Users/femal/Dropbox/Soil_C_modeling/Modern/Data")

#extract relevant values at sites

raw_data = read_xlsx("Modern_Data_Final.xlsx", sheet= 1)
sites = read_xlsx("Modern_Data_Final.xlsx", sheet = 2)
coords = matrix(c(sites$Lon, sites$Lat), nrow(sites), 2)

sites$Tma = extract(Tma, coords)
sites$Pa = extract(Pa, coords)
sites$TmWQ_min_a = extract(TmWQ_min_a, coords)
sites$PfWQ = extract(PfWQ, coords)
sites$TmDQ_min_a = extract(TmDQ_min_a, coords)
sites$PfDQ = extract(PfDQ, coords)
sites$DQ = extract(DQ, coords)
sites$Ra = extract(Ra, coords)
sites$per.C4 = extract(C4, coords)
sites[,1] <- NULL
write.csv(sites, "Modern_Sites_Final.csv")

# Isotope data averaged over all depths

data.aves = aggregate.data.frame(raw_data, list(raw_data$Site), mean)
data.aves$Site <- data.aves$Group.1
data.aves[,1] <- NULL
data.aves[,1] <- NULL
data.aves[,2] <- NULL
View(data.aves)
write.csv(data.aves, "Modern_data_aves_Final.csv")

sites.depth <- merge.data.frame(raw_data, sites, by.x="Site", by.y="Site", all.x=T)

# Delete any sites that don't have depth info

sites.depth <- subset(sites.depth, !is.na(sites.depth$depth.cm))
length(unique(sites.depth$Site))
write.csv(sites.depth, "Modern_Sites_depth_Final.csv")

# If you already have done this, read in data.aves and sites

setwd("C:/Users/femal/Dropbox/Soil_C_modeling/Modern/Data/")

data.aves = read.csv("Modern_data_aves_Final.csv")
sites = read.csv("Modern_Sites_Final.csv")
sites.depth = read.csv("Modern_Sites_depth_Final.csv")
raw_data = read_xlsx("Modern_Data_Final.xlsx", sheet= 1)

# Site isotope sds
raw_data = as.data.frame(raw_data)

site.sds = data.frame(Site = character(0), d13C_sd = numeric(0), d18O_sd = numeric(0))
for(i in unique(raw_data$Site)){
  
  if(length(raw_data[which(raw_data$Site == i),"d13C.measured"]) > 2){
    
    site.sd = data.frame("Site" = i, "d13C_sd" = sd(raw_data[which(raw_data$Site == i),"d13C.measured"]), "d18O_sd" = sd(raw_data[which(raw_data$Site == i),"d18O.measured"]))
    site.sds = rbind(site.sds, site.sd)
    
  }}

O_site_sd = mean(na.omit(site.sds$d18O_sd))
C_site_sd = mean(na.omit(site.sds$d13C_sd))
O_site_sd
C_site_sd

View(sites)

setwd("C:/Users/femal/Dropbox/Soil_C_modeling/Modern/R Scripts/")

#Run theoretical model
source("Theoretical.R")

#Run theoretical model with depth as input
source("Theoretical_depth.R")

#Run Seasonal Precip and Evap Optimization
source("Evap_Seasonal_OPT.R")

#Run Seasonal Precip and Evap Optimization with depth as input
source("Evap_Seasonal_OPT_depth.R")

#Run respiration optimization
source("Respiration_OPT.R")

#Run respiration optmiziation with depth as input
source("Respiration_OPT_depth.R")

#Apply optimization
source("Optimized.R")

#Apply optimization with depth as input
source("Optimized_depth.R")

#Sensitivity tests
source("SensTest.R")

#Prediction maps - worldwide
source("PredMaps_OPT.R")

setwd("C:/Users/femal/Dropbox/Soil_C_modeling/Modern/R Scripts/")

save.image("2020_Final.RData")
