library(raster)
library(rgdal)

setwd("C:/Users/femal/Dropbox/worldclim/")

TmQ_djf <- raster("TmQ_djf.tif")
TmQ_mam <- raster("TmQ_mam.tif")
TmQ_jja <- raster("TmQ_jja.tif")
TmQ_son <- raster("TmQ_son.tif")

setwd("C:/Users/femal/Dropbox/Soil_C_modeling/Modern/Data")
d18O_Tma = read.csv("Oma.csv")
sites_OIPC <- read.csv("Modern_Sites_OIPC_Final.csv")
coords = matrix(c(sites_OIPC $Lon, sites_OIPC $Lat), nrow(sites_OIPC), 2)
sites_OIPC$Tma = extract(Tma, coords)
sites_OIPC$Pa = extract(Pa, coords)
sites_OIPC$TmWQ_min_a = extract(TmWQ_min_a, coords)
sites_OIPC$PfWQ = extract(PfWQ, coords)
sites_OIPC$TmDQ_min_a = extract(TmDQ_min_a, coords)
sites_OIPC$PfDQ = extract(PfDQ, coords)
sites_OIPC$DQ = extract(DQ, coords)
sites_OIPC$Ra = extract(Ra, coords)
sites_OIPC$per.C4 = extract(C4, coords)
sites_OIPC$TmQ_djf = extract(TmQ_djf, coords)
sites_OIPC$TmQ_mam = extract(TmQ_mam, coords)
sites_OIPC$TmQ_jja = extract(TmQ_jja, coords)
sites_OIPC$TmQ_son = extract(TmQ_son, coords)
coords = matrix(c(d18O_Tma$Longitude, d18O_Tma$Latitude), nrow(d18O_Tma), 2)
d18O_Tma$Tma = extract(Tma, coords)

sites_OIPC$d18O_OIPC_djf = rowMeans(cbind(sites_OIPC[,30], sites_OIPC[,19], sites_OIPC[,20]))
sites_OIPC$d18O_OIPC_mam = rowMeans(sites_OIPC[,21:23])
sites_OIPC$d18O_OIPC_jja = rowMeans(sites_OIPC[,24:26])
sites_OIPC$d18O_OIPC_son = rowMeans(sites_OIPC[,26:29])

d18O_lm <- lm(d18O_Tma$d18O ~ d18O_Tma$Tma)
summary(d18O_lm)
d18O_Tma_20_80 = subset(d18O_Tma, abs(d18O_Tma$Latitude) > 20 & abs(d18O_Tma$Latitude) < 80 & d18O_Tma$nyrs > 1)
d18O_lm_20_80 = lm(d18O_Tma_20_80$d18O ~ d18O_Tma_20_80$Tma)
summary(d18O_lm_20_80)
err <- d18O_Tma_20_80$d18O - (d18O_Tma_20_80$Tma * 0.60 - 15.57)
err = na.omit(err)
length(err)
d18O_rmse <- sqrt(mean(err^2)) 
d18O_rmse

sites_OIPC$d18O_Model_ann <- sites_OIPC$Tma * 0.6 - 15.57
sites_OIPC$d18O_Model_djf <- sites_OIPC$TmQ_djf * 0.6 - 15.57
sites_OIPC$d18O_Model_mam <- sites_OIPC$TmQ_mam * 0.6 - 15.57
sites_OIPC$d18O_Model_jja <- sites_OIPC$TmQ_jja * 0.6 - 15.57
sites_OIPC$d18O_Model_son <- sites_OIPC$TmQ_son * 0.6 - 15.57

plot(sites_OIPC$d18O_OIPC_djf ~ sites_OIPC$d18O_Model_djf)

d18O_OIPC_allseas <- c(sites_OIPC$d18O_OIPC_djf, sites_OIPC$d18O_OIPC_mam, sites_OIPC$d18O_OIPC_jja, sites_OIPC$d18O_OIPC_son)
d18O_Model_allseas <- c(sites_OIPC$d18O_Model_djf, sites_OIPC$d18O_Model_mam, sites_OIPC$d18O_Model_jja, sites_OIPC$d18O_Model_son)

allseas_lm = lm(d18O_OIPC_allseas ~ d18O_Model_allseas)
summary(allseas_lm)
err = d18O_OIPC_allseas - d18O_Model_allseas
allseas_rmse <- sqrt(mean(err^2)) 
allseas_rmse
meanann_lm = lm(sites_OIPC$d18O_Model_ann ~ sites_OIPC$d18O_OIPC_ann)
summary(meanann_lm)
err = sites_OIPC$d18O_Model_ann - sites_OIPC$d18O_OIPC_ann
meanann_rmse <- sqrt(mean(err^2)) 
meanann_rmse

