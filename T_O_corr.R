library(raster)
library(rgdal)

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

