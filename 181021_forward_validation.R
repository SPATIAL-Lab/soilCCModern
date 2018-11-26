###Code for conducting validation of soil carbonate model

setwd("~/GitHub/soilCCModern")

#Prepare gridded climate data
library(raster)

#read in worldclim 2.5min mean temperature grids, downloaded 9-18-18
setwd("C:/Users/gjbowen/Dropbox/Archived/ArcGIS/worldclim/")
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

#####That's all the grid processing - next bit needs only be run once after data are updated

#read in the compiled data
setwd("C:/Users/gjbowen/Dropbox/HypoMirror/soilCCModern/")
library(xlsx)

#extract relevant values at sites
precipcomp <- read.csv("valsites_sel_d18O_P.csv")
rawsites = read.xlsx("modern_comparison.xlsx", sheetIndex = 3)
sites = read.xlsx("modern_comparison_raw.xlsx", sheetIndex = 1)
coords = matrix(c(sites$Lon, sites$Lat), nrow(sites), 2)
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

######Now code to run some simulations
nsynth = 5000

sites = read.csv("valsites_sel.csv")

#read in real data
data = read.xlsx("modern_comparison.xlsx", sheetIndex = 3)
data.aves = aggregate.data.frame(data, list(data$Site), mean, simplify = TRUE)
data.aves$Site = NULL
data.comp = merge.data.frame(sites, data.aves, by.x = "Site", by.y = "Group.1")

## Compare dO_P from OIPC and model
precipcomp <- merge.data.frame(precipcomp, hq.comp, by.x="Site", by.y="Site")

c = ceiling((precipcomp$map.wc.y / max(precipcomp$map.wc.y)) * 5)

pal = rainbow(5)

layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow=T))

plot(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_ann, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC Ann",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Model Ann",delta^{18}, "O (\u2030)")))
points(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_ann, pch=1)
abline(0,1)
legend("bottomright", title = "MAP (mm)", fill=rainbow(5), legend=c("0 - 147", "147 - 294","294 - 441", "441 - 588" , "588 - 737"))

plot(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_jja, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC JJA",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Model Ann",delta^{18}, "O (\u2030)")))
points(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_jja, pch=1)
abline(0,1)

plot(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_mam, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC MAM",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Model Ann",delta^{18}, "O (\u2030)")))
points(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_mam, pch=1)
abline(0,1)

plot(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_son, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC SON",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Model Ann",delta^{18}, "O (\u2030)")))
points(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_son, pch=1)
abline(0,1)

plot(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_djf, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC DJF",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Model Ann",delta^{18}, "O (\u2030)")))
points(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_djf, pch=1)
abline(0,1)


plot(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_son, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC SON",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Model Ann",delta^{18}, "O (\u2030)")))
points(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_son, pch=1)
abline(0,1)


### Compare measured d18O carbonate (w/ no evap) in HQ and DQ temps to OIPC d18O precip

plot(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_djf, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC DJF ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure HQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_djf, pch=1)
abline(0,1)
legend("bottomright", title = "MAP (mm)", fill=rainbow(5), legend=c("0 - 147", "147 - 294","294 - 441", "441 - 588" , "588 - 737"))


plot(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_jja, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC JJA ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure HQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_jja, pch=1)
abline(0,1)

plot(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_mam, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC MAM ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure HQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_mam, pch=1)
abline(0,1)

plot(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_son, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC SON ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure HQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_son, pch=1)
abline(0,1)

plot(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_ann, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC Ann ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure HQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_ann, pch=1)
abline(0,1)

plot(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_Model_ann, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("Model Ann ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure HQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_Model_ann, pch=1)
abline(0,1)


## DQ
plot(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_djf, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC DJF ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure DQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_djf, pch=1)
abline(0,1)
legend("bottomright", title = "MAP (mm)", fill=rainbow(5), legend=c("0 - 147", "147 - 294","294 - 441", "441 - 588" , "588 - 737"))


plot(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_jja, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC JJA ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure DQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_jja, pch=1)
abline(0,1)

plot(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_mam, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC MAM ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure DQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_mam, pch=1)
abline(0,1)

plot(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_son, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC SON ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure DQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_son, pch=1)
abline(0,1)

plot(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_ann, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC Ann ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure DQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_ann, pch=1)
abline(0,1)

plot(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_Model_ann, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("Model Ann ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure DQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_Model_ann, pch=1)
abline(0,1)

## HQ Vs DQ dO_P
plot(precipcomp$dO_P_carb_dq ~ precipcomp$dO_P_carb_hq, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("Measure HQ ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure DQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_dq ~ precipcomp$dO_P_carb_hq, pch=1)
abline(0,1)
legend("bottomright", title = "MAP (mm)", fill=rainbow(5), legend=c("0 - 147", "147 - 294","294 - 441", "441 - 588" , "588 - 737"))


c = ceiling((precipcomp$Alt.y / max(precipcomp$Alt.y)) * 5)

pal = rainbow(5)

layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow=T))

plot(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_ann, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC Ann",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Model Ann",delta^{18}, "O (\u2030)")))
points(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_ann, pch=1)
abline(0,1)
legend("bottomright", title = "Alt (m)", fill=rainbow(5), legend=c("0 - 500", "500 - 1000","1000 - 1500", "1500 - 2000" , "2000 - 2500"))

plot(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_jja, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC JJA",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Model Ann",delta^{18}, "O (\u2030)")))
points(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_jja, pch=1)
abline(0,1)

plot(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_mam, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC MAM",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Model Ann",delta^{18}, "O (\u2030)")))
points(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_mam, pch=1)
abline(0,1)

plot(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_son, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC SON",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Model Ann",delta^{18}, "O (\u2030)")))
points(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_son, pch=1)
abline(0,1)

plot(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_djf, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC DJF",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Model Ann",delta^{18}, "O (\u2030)")))
points(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_djf, pch=1)
abline(0,1)


plot(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_son, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC SON",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Model Ann",delta^{18}, "O (\u2030)")))
points(precipcomp$d18O_Model_ann ~ precipcomp$d18O_OIPC_son, pch=1)
abline(0,1)


### Compare measured d18O carbonate (w/ no evap) in HQ and DQ temps to OIPC d18O precip

plot(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_djf, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC DJF ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure HQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_djf, pch=1)
abline(0,1)
legend("bottomright", title = "Alt (m)", fill=rainbow(5), legend=c("0 - 500", "500 - 1000","1000 - 1500", "1500 - 2000" , "2000 - 2500"))


plot(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_jja, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC JJA ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure HQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_jja, pch=1)
abline(0,1)

plot(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_mam, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC MAM ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure HQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_mam, pch=1)
abline(0,1)

plot(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_son, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC SON ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure HQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_son, pch=1)
abline(0,1)

plot(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_ann, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC Ann ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure HQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_OIPC_ann, pch=1)
abline(0,1)

plot(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_Model_ann, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("Model Ann ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure HQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_hq ~ precipcomp$d18O_Model_ann, pch=1)
abline(0,1)


## DQ
plot(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_djf, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC DJF ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure DQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_djf, pch=1)
abline(0,1)
legend("bottomright", title = "Alt (m)", fill=rainbow(5), legend=c("0 - 500", "500 - 1000","1000 - 1500", "1500 - 2000" , "2000 - 2500"))


plot(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_jja, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC JJA ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure DQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_jja, pch=1)
abline(0,1)

plot(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_mam, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC MAM ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure DQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_mam, pch=1)
abline(0,1)

plot(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_son, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC SON ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure DQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_son, pch=1)
abline(0,1)

plot(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_ann, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("OIPC Ann ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure DQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_OIPC_ann, pch=1)
abline(0,1)

plot(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_Model_ann, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("Model Ann ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure DQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_dq ~ precipcomp$d18O_Model_ann, pch=1)
abline(0,1)

## HQ Vs DQ dO_P
plot(precipcomp$dO_P_carb_dq ~ precipcomp$dO_P_carb_hq, pch=16, col=pal[c], xlim=c(-20,0), ylim=c(-20,0),
     xlab=expression(paste("Measure HQ ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Measure DQ ",delta^{18}, "O (\u2030)")))
points(precipcomp$dO_P_carb_dq ~ precipcomp$dO_P_carb_hq, pch=1)
abline(0,1)
legend("bottomright", title = "Alt (m)", fill=rainbow(5), legend=c("0 - 500", "500 - 1000","1000 - 1500", "1500 - 2000" , "2000 - 2500"))


dev.off()
## Estimate oxygen isotopes of precip from d18O carbonate and hq and dq temperatures

for (i in 1:nrow(precipcomp)){
  
A_O <- 2.71828 ^ ((2.78e6 / (precipcomp[i,"mat.wc.y"] + precipcomp[i, "hqt.offset.y"] + 273.15) ^ 2 - 2.89) / 1000)
R_O <- (precipcomp[i, "d18O.measured"] / 1000 + 1) * 0.002067181
R_O_P <-R_O / A_O 
precipcomp[i, "dO_P_carb_hq"] <-  (R_O_P / 0.0020052 - 1) * 1000

A_O <- 2.71828 ^ ((2.78e6 / (precipcomp[i, "mat.wc.y"] + precipcomp[i, "dqt.offset.y"] + 273.15) ^ 2 - 2.89) / 1000)
R_O <- (precipcomp[i, "d18O.measured"] / 1000 + 1) * 0.002067181
R_O_P <- R_O / A_O 
precipcomp[i, "dO_P_carb_dq"] <-  (R_O_P / 0.0020052 - 1) * 1000

}

### This part runs only sites w/ clumpted data
data.clump = data.comp[!is.na(data.comp$D47.measured),]
for(i in 1:nrow(data.clump)){
  if(is.na(data.clump$MAT[i])){data.clump$MAT[i] = data.clump$mat.wc[i]}
  if(is.na(data.clump$MAP[i])){data.clump$MAP[i] = data.clump$map.wc[i]}
}

clump_pred_meantemp = data.frame(depth=numeric(0), soil18O=numeric(0), d13C=numeric(0), d18O=numeric(0))

## Calibration curve of Defliese et. al. 2017. 90 degree rxn only
data.clump$Clumped_T = sqrt(0.0422e6 / (data.clump$D47.measured - 0.1262)) - 273.15 

## Calibration curve of Kelson et. al. 2017
data.clump$Clumped_T = sqrt(0.0422e6 / (data.clump$D47.measured - 0.215)) - 273.15 

## Calibration curve of Ghosh 2006, recalibrated by Dennis 2011
data.clump$Clumped_T = sqrt(0.0636e6 / (data.clump$D47.measured + 0.0047)) - 273.15 

## Calibration curve of Dennis and Schrag 2010, recalibrated by Dennis 2011
data.clump$Clumped_T = sqrt(0.0362e6 / (data.clump$D47.measured - 0.292)) - 273.15 

## Is the clumped temperature closer to hqt or dqt?
data.clump$Clumped.offset = data.clump$Clumped_T - data.clump$mat.wc
data.clump$Clumped.offsetmean = (data.clump$hqt.offset + data.clump$dqt.offset) / 2


for(i in 1:nrow(data.clump)){
  
if(data.clump$Clumped.diff[i] > data.clump$Clumped.offsetmean[i]) {
                              data.clump$Clumped.frac[i] = data.clump$hqp.frac[i] 
                              data.clump$Informed.offset[i] = data.clump$hqt.offset[i]} 
                              else {
                              data.clump$Clumped.frac[i] = data.clump$dqp.frac[i] 
                              data.clump$Informed.offset[i] = data.clump$dqt.offset[i]}
  }

data.clump$Informed.Temp = data.clump$Informed.offset + data.clump$mat.wc

## Running it in the forward model - mean temperature in the "informed by clumped" quarter
for(i in 1:nrow(data.clump)){
  clump_pred_meantemp[i,] = sm_forward(data.clump$map.wc[i], data.clump$mat.wc[i], data.clump$Clumped.frac[i], data.clump$Informed.offset[i], 280)
  
}
clump_pred_meantemp$Site = data.clump$Site

# Add it to the predictions and plot

clump.comp.meantemp = merge.data.frame(clump_pred_meantemp, data.clump, by.x = "Site", by.y = "Site", all.x=TRUE)

c = ceiling((clump.comp.meantemp$map.wc / max(clump.comp.meantemp$map.wc)) * 5)

pal = rainbow(5)

layout(matrix(c(1,2), 1, 2, byrow=T))

plot(clump.comp.meantemp$d13C, clump.comp.meantemp$d13C.measured, pch=16, col=pal[c], xlim=c(-14,0), ylim=c(-14,0), main="Using Informed Mean Temperature",
     xlab=expression(paste("Predicted ",delta^{13}, "C (\u2030)")),
     ylab=expression(paste("Observed ",delta^{13}, "C (\u2030)")))
points(clump.comp.meantemp$d13C, clump.comp.meantemp$d13C.measured, pch=1)
abline(0,1)
legend("bottomright", title = "MAP (mm)", fill=rainbow(5), legend=c("0 - 136", "136 - 272", "272 - 408" , "408 - 544", "544 - 680"))


plot(clump.comp.meantemp$d18O, clump.comp.meantemp$d18O.measured, pch=16, col=pal[c], xlim=c(-16,0), ylim=c(-16,0),
     xlab=expression(paste("Predicted ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Observed ",delta^{18}, "O (\u2030)")))
points(clump.comp.meantemp$d18O, clump.comp.meantemp$d18O.measured, pch=1)
abline(0,1)

## Running it in the forward model - using clumped temperature as the mean quarter temperature
clump_pred = data.frame(depth=numeric(0), soil18O=numeric(0), d13C=numeric(0), d18O=numeric(0))

for(i in 1:nrow(data.clump)){
  clump_pred[i,] = sm_forward(data.clump$map.wc[i], data.clump$mat.wc[i], data.clump$Clumped.frac[i], data.clump$Clumped.diff[i], 280)
  
}
clump_pred$Site = data.clump$Site

# Add it to the predictions and plot
clump.comp = merge.data.frame(clump_pred, data.clump, by.x = "Site", by.y = "Site", all.x=TRUE)

layout(matrix(c(1,2), 1, 2, byrow=T))
plot(clump.comp$d13C, clump.comp$d13C.measured, xlim=c(-14,0), ylim=c(-14,0), main="Using Clumped Temperature Directly")
abline(0,1)

plot(clump.comp$d18O, clump.comp$d18O.measured, xlim=c(-16,-1), ylim=c(-16,-1))
abline(0,1)


# Temperature comparison Graphs

layout(matrix(c(1,2,3), 1, 3, byrow=T))
plot(clump.comp.meantemp$Informed.Temp, clump.comp.meantemp$Clumped_T, ylim = c(0,56), xlim=c(0,56), xlab="Informed Mean Temperature (C)", ylab="Clumped Temperature (C)")
abline(0,1)

plot(clump.comp.meantemp$mat.wc + clump.comp.meantemp$hqt.offset, clump.comp.meantemp$Clumped_T, ylim = c(0,56), xlim=c(0,56), xlab="Hot Quarter Temperature (C)", ylab="Clumped Temperature (C)")
abline(0,1)

plot(clump.comp.meantemp$mat.wc + clump.comp.meantemp$dqt.offset, clump.comp.meantemp$Clumped_T, ylim = c(0,56), xlim=c(0,56), xlab="Dry Quarter Temperature (C)", ylab="Clumped Temperature (C)")
abline(0,1)

# Temperature comparison stats
mean(clump.comp.meantemp$hqt.offset)
sd(clump.comp.meantemp$hqt.offset)
mean(clump.comp.meantemp$dqt.offset)
sd(clump.comp.meantemp$dqt.offset)
mean(clump.comp.meantemp$Clumped.offset)
sd(clump.comp.meantemp$Clumped.offset)
mean(clump.comp.meantemp$hqt.offset - clump.comp.meantemp$Clumped.offset)
sd(clump.comp.meantemp$hqt.offset - clump.comp.meantemp$Clumped.offset)
mean(clump.comp.meantemp$dqt.offset - clump.comp.meantemp$Clumped.offset)
sd(clump.comp.meantemp$dqt.offset - clump.comp.meantemp$Clumped.offset)

## Exclude hyperarid sites
clump.comp.sel = subset(clump.comp, clump.comp$MAP > 100)
plot(clump.comp.sel$d13C, clump.comp.sel$d13C.measured)
abline(0,1)
plot(clump.comp.sel$d18O, clump.comp.sel$d18O.measured)
abline(0,1)

## Respiration data subset into subhumid to arid

srdb = read.csv("srdb-data.csv")
srdb_carb <- subset(srdb, srdb$MAP < 300)

## Subset of sites that have depth information - matching them with the depth that is closest (w/in 10 cm) of the model estimated depths
depth = data.frame(site=numeric(0), depthM = numeric(0), depthObs = numeric(0), d13CD = numeric(0), d18OD = numeric(0))
for(i in 1:nrow(hq.comp)) {
  for(j in 1:nrow(rawsites)){
if(identical(as.character(hq.comp[i,1]), as.character(rawsites[j,2])) && is.na(rawsites[j,3]) == F && hq.comp[i,2] > rawsites[j,3] - 10 && hq.comp[i,2] < rawsites[j,3] + 10) 
{depthCO = data.frame(site=rawsites[j,2], depthM = hq.comp[i,2], depthObs = rawsites[j,3], d13C_measured=rawsites[j,4], d18O_measured=rawsites[j,5])
depth = rbind(depth, depthCO)}
  }}
## Average each site
depth.aves = aggregate.data.frame(depth, list(depth$site), mean, simplify = TRUE)
depth.aves$site=NULL
# Add other data
sites.depth <- merge.data.frame(depth.aves, sites, by.x = "Group.1", by.y = "Site")

## Run fw model for only depth sites - hq and dq w. evap

## w/ evap
hq_pred.depth = data.frame(depth=numeric(0), soil18O=numeric(0), dO_P=numeric(0), d13C=numeric(0), d18O=numeric(0))
for(i in 1: nrow(sites.depth)){
  hq_pred.depth[i,] = sm_forward_evap(sites.depth$map.wc[i], sites.depth$mat.wc[i], sites.depth$hqp.frac[i], sites.depth$hqt.offset[i], 280)
  
}
hq_pred.depth$Site = sites.depth$Group.1

dq_pred.depth = data.frame(depth=numeric(0), soil18O=numeric(0), dO_P=numeric(0),d13C=numeric(0), d18O=numeric(0))
for(i in 1: nrow(sites.depth)){
  dq_pred.depth[i,] = sm_forward_evap(sites.depth$map.wc[i], sites.depth$mat.wc[i], sites.depth$dqp.frac[i], sites.depth$dqt.offset[i], 280)
  
}
dq_pred.depth$Site = sites.depth$Group.1


layout(matrix(c(1,2,3,4), 2, 2, byrow=T))

hq.comp.depth = merge.data.frame(hq_pred.depth, sites.depth, by.x = "Site", by.y = "Group.1", all.x=TRUE)
plot(hq.comp.depth$d13C, hq.comp.depth$d13C_measured)
abline(0,1)
plot(hq.comp.depth$d18O, hq.comp.depth$d18O_measured)
abline(0,1)

dq.comp.depth = merge.data.frame(dq_pred.depth, sites.depth, by.x = "Site", by.y = "Group.1")
plot(dq.comp.depth$d13C, dq.comp.depth$d13C_measured)
abline(0,1)
plot(dq.comp.depth$d18O, dq.comp.depth$d18O_measured)
abline(0,1)

# Same results as the overall averages w/o depth

##
#run this for only the not superarid ones
sites = read.csv("valsites_sel.csv")


sites = sites[sites$map.wc > 100,]

hq_pred = data.frame(depth=numeric(0), soil18O=numeric(0), d13C=numeric(0), d18O=numeric(0))
for(i in 1: nrow(sites)){
  hq_pred[i,] = sm_forward(sites$map.wc[i], sites$mat.wc[i], sites$hqp.frac[i], sites$hqt.offset[i], 280)
  
}
hq_pred$Site = sites$Site

dq_pred = data.frame(depth=numeric(0), soil18O=numeric(0), d13C=numeric(0), d18O=numeric(0))
for(i in 1: nrow(sites)){
  dq_pred[i,] = sm_forward(sites$map.wc[i], sites$mat.wc[i], sites$dqp.frac[i], sites$dqt.offset[i], 280)
  
}
dq_pred$Site = sites$Site

## w/ evap
hq_pred = data.frame(depth=numeric(0), soil18O=numeric(0), dO_P=numeric(0), d13C=numeric(0), d18O=numeric(0))
for(i in 1: nrow(sites)){
  hq_pred[i,] = sm_forward_evap(sites$map.wc[i], sites$mat.wc[i], sites$hqp.frac[i], sites$hqt.offset[i], 280)
  
}
hq_pred$Site = sites$Site

dq_pred = data.frame(depth=numeric(0), soil18O=numeric(0), d13C=numeric(0), d18O=numeric(0))
for(i in 1: nrow(sites)){
  dq_pred[i,] = sm_forward_evap(sites$map.wc[i], sites$mat.wc[i], sites$dqp.frac[i], sites$dqt.offset[i], 280)
  
}
dq_pred$Site = sites$Site

#add it to the predictions and plot
layout(matrix(c(1,2,3,4), 2, 2, byrow=T))

hq.comp = merge.data.frame(hq_pred, data.aves, by.x = "Site", by.y = "Group.1", all.x=TRUE)
plot(hq.comp$d13C, hq.comp$d13C.measured)
abline(0,1)
plot(hq.comp$d18O, hq.comp$d18O.measured)
abline(0,1)

dq.comp = merge.data.frame(dq_pred, data.aves, by.x = "Site", by.y = "Group.1")
plot(dq.comp$d13C, dq.comp$d13C.measured)
abline(0,1)
plot(dq.comp$d18O, dq.comp$d18O.measured)
abline(0,1)

## Plots w/ MAP colors
hq.comp = merge.data.frame(hq.comp, sites, by.x = "Site", by.y = "Site", all.x=TRUE)
dq.comp = merge.data.frame(dq.comp, sites, by.x = "Site", by.y = "Site", all.x=TRUE)

c = ceiling((hq.comp$map.wc / max(hq.comp$map.wc)) * 6)

pal = rainbow(6)

layout(matrix(c(1,2,3,4), 2, 2, byrow=T))

plot(hq.comp$d13C, hq.comp$d13C.measured, pch=16, col=pal[c], xlim=c(-14,0), ylim=c(-14,0), main="HQ Carbon", cex = 1.25,
     xlab=expression(paste("Predicted ",delta^{13}, "C (\u2030)")),
     ylab=expression(paste("Observed ",delta^{13}, "C (\u2030)")))
points(hq.comp$d13C, hq.comp$d13C.measured, pch=1)
abline(0,1)
legend("bottomright", title = "MAP (mm)", fill=rainbow(6), legend=c("0-122", "122-244", "244-366", "366-488", "488-600", "600+"))


plot(hq.comp$d18O, hq.comp$d18O.measured, pch=16, col=pal[c], xlim=c(-16,0), ylim=c(-16,0),main="HQ Oxygen w/ Evap",cex = 1.25,
     xlab=expression(paste("Predicted ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Observed ",delta^{18}, "O (\u2030)")))
points(hq.comp$d18O, hq.comp$d18O.measured, pch=1)
abline(0,1)

plot(dq.comp$d13C, dq.comp$d13C.measured, pch=16, col=pal[c], xlim=c(-14,0), ylim=c(-14,0), main="DQ Carbon",cex = 1.25,
     xlab=expression(paste("Predicted ",delta^{13}, "C (\u2030)")),
     ylab=expression(paste("Observed ",delta^{13}, "C (\u2030)")))
points(dq.comp$d13C, dq.comp$d13C.measured, pch=1)
abline(0,1)

plot(dq.comp$d18O, dq.comp$d18O.measured, pch=16, col=pal[c], xlim=c(-16,0), ylim=c(-16,0), main="DQ Oxygen w/ Evap",cex = 1.25,
     xlab=expression(paste("Predicted ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Observed ",delta^{18}, "O (\u2030)")))
points(dq.comp$d18O, dq.comp$d18O.measured, pch=1)
abline(0,1)

## Calculate RMSE

  C.residuals <- hq.comp$d13C.measured - hq.comp$d13C
  O.residuals <- hq.comp$d18O.measured - hq.comp$d18O
  
  C.rmse <- sqrt(mean(abs(C.residuals^2)))
  O.rmse <- sqrt(mean(abs(O.residuals^2)))

  C.rmse
  O.rmse
  
############### Forward model function for use in sensitivity testing

sm_forward = function(MAP, MAT, P_seas, T_seas, pCO2){

    deltaA = rnorm(nsynth, -6.5, 0.3)
    pores = rnorm(nsynth, 0.46, 0.1)
    tort = rnorm(nsynth, 0.7, 0.1)
    
    #Solar radiation, here fixed
    Rs = 20.35
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

    #Convert pCO2 to units mol/cm^3
    pCO2_mcc = pCO2 / (0.08206 * CQT_K * 10^9)  #mol/cm^3

    #Relative humidity, now beta distribution
    h_m <- 0.25 + 0.7 * (CQP / 900)
    h_var = 0.05^2
    size = h_m*(1-h_m)/h_var - 1
    alpha = h_m * size
    beta = (1-h_m) * size
    h = rbeta(nsynth, alpha, beta)
    RH <- h * 100
    
    # Precipitation O isotope ratios 
    dO_P_m <- -13.7 + 0.55 * MAT
    dO_P = rnorm(nsynth, dO_P_m, 1.7)
    R_O_P = (dO_P / 1000 + 1) * RO.vsmow

    #Depth to carbonate, now gamma dist
    z_mean <- MAP*0.0925 + 13.4  #top of Bk equation based on Retallack 2005 data
    z_thick = abs(CMP_mm - MAP/12) * 0.74 + 17.3 #thickness of Bk, CMP as proxy for seasonality 
    z_mean = z_mean + z_thick/2 #find middle of Bk in cm
    theta = 20^2/z_mean  #gamma scale parameter, using 20cm as variance from Retallack fit
    k = z_mean / theta #gamma shape parameter
    z = rgamma(nsynth, shape = k, scale = theta)
    z <- pmin(z, 100) 
    
    #Respiration rate, now gamma dist
    R_day_m <- 1.25 * exp(0.0545 * CQT) * CMP_cm / (4.259 + CMP_cm)  #Raich 2002, gC/m2day
    R_day_m <- 1.24 * exp(0.055 * CQT) * CMP_cm / (4.78 + CMP_cm) # Re-parameterized to MAP < 760, switch equations for comparison
    theta = (R_day_m*0.5)^2/R_day_m #gamma scale parameter, using mean residual of 50% based on Raich validation data 
    k = R_day_m / theta #gamma shape parameter
    R_day = rgamma(nsynth, shape = k, scale = theta) #lets use gamma for these quants bounded at zero....
    R_day = R_day / (12.01 * 100^2)  #molC/cm2day
    R_sec <- R_day / (24 * 3600)  #molC/cm2s
    R_sec = R_sec / L / pores #molC/cm3s
    
    #
    
    #Potential ET
    ETP_D_m <- ifelse (RH < 50, 0.0133 * (CQT / (CQT + 15)) * (1/23.885 * Rs + 50) * (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (1/23.885 * Rs + 50))
    ETP_D = rnorm(nsynth, ETP_D_m, 0.2)  #PET in mm/day, Turc 1961
    ETP_M <- ETP_D * 30  #mm/month
    
    #Actual ET
    ETA_var = rnorm(nsynth, 1, 0.2) #This noise parmeter limits ETA<CMP_mm but allows variation around ETP, as observed
    ETA = CMP_mm*3 * (1 / (sqrt(1 + (1 / ((ETP_M / (CMP_mm*3)) * ETA_var)) ^ 2))) #AET in mm/month from Budyko curve
    #here scaled eta to quarter precip, assuming potential carry-over
    
    #Free air porosity
    #Have updated, now scales volumetrically w/ excess precipitation relative to pore space
    FAP <- pmin((pores - (CMP_mm - ETA)/(L*10*pores)), pores)
    FAP = pmax(FAP,0.01) #dimensionless
    
    #CO2 Diffusion coefficients
    DIFC = FAP * tort * 0.1369 * (CQT_K / 273.15) ^ 1.958

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
    
    #Soil carbonate C isotopes
    A_CO2_Carb <- 2.71828 ^ (-2.988e3 / CQT_K ^ 2 + 7.6663 / CQT_K - 0.0024612)
    R_Soil <- (dC_Soil / 1000 + 1) * RC.vpdb
    R_Carb <- R_Soil / A_CO2_Carb

    #Soil carbonate O isotopes
    A_O <- 2.71828 ^ ((2.78e6 / CQT_K ^ 2 - 2.89) / 1000)
    R_O_Carb <- R_O_P * A_O
    
    dC_Carb <- (R_Carb / RC.vpdb - 1) * 1000 
    dO_Carb <- (R_O_Carb / RO.vpdb - 1) * 1000
  
    dat = c(median(dC_Carb), median(dO_Carb))
  
  return(dat)
}


sm_forward_evap = function(MAP, MAT, P_seas, T_seas, pCO2){
  
  deltaA = rnorm(nsynth, -6.5, 0.3)
  pores = rnorm(nsynth, 0.46, 0.1)
  tort = rnorm(nsynth, 0.7, 0.1)
  esw = 1
  
  #Solar radiation, here fixed
  Rs = 20.35
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
  
  #Convert pCO2 to units mol/cm^3
  pCO2_mcc = pCO2 / (0.08206 * CQT_K * 10^9)  #mol/cm^3
  
  #Relative humidity, now beta distribution
  h_m <- 0.25 + 0.7 * (CQP / 900)
  h_var = 0.05^2
  size = h_m*(1-h_m)/h_var - 1
  alpha = h_m * size
  beta = (1-h_m) * size
  h = rbeta(nsynth, alpha, beta)
  RH <- h * 100
  
  # Precipitation O isotope ratios 
  dO_P_m <- -13.7 + 0.55 * MAT
  dO_P = rnorm(nsynth, dO_P_m, 1.7)
  R_O_P = (dO_P / 1000 + 1) * RO.vsmow
  
  #Depth to carbonate, now gamma dist
  z_mean <- MAP*0.0925 + 13.4  #top of Bk equation based on Retallack 2005 data
  z_thick = abs(CMP_mm - MAP/12) * 0.74 + 17.3 #thickness of Bk, CMP as proxy for seasonality 
  z_mean = z_mean + z_thick/2 #find middle of Bk in cm
  theta = 20^2/z_mean  #gamma scale parameter, using 20cm as variance from Retallack fit
  k = z_mean / theta #gamma shape parameter
  z = rgamma(nsynth, shape = k, scale = theta)
  z <- pmin(z, 100) 
  
  #Respiration rate, now gamma dist
  R_day_m <- 1.25 * exp(0.0545 * CQT) * CMP_cm / (4.259 + CMP_cm)  #Raich 2002, gC/m2day
  R_day_m <- 1.24 * exp(0.055 * CQT) * CMP_cm / (4.78 + CMP_cm) # Re-parameterized to MAP < 760, switch equations for comparison
  theta = (R_day_m*0.5)^2/R_day_m #gamma scale parameter, using mean residual of 50% based on Raich validation data 
  k = R_day_m / theta #gamma shape parameter
  R_day = rgamma(nsynth, shape = k, scale = theta) #lets use gamma for these quants bounded at zero....
  R_day = R_day / (12.01 * 100^2)  #molC/cm2day
  R_sec <- R_day / (24 * 3600)  #molC/cm2s
  R_sec = R_sec / L / pores #molC/cm3s
  
  #Potential ET
  ETP_D_m <- ifelse (RH < 50, 0.0133 * (CQT / (CQT + 15)) * (1/23.885 * Rs + 50) * (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (1/23.885 * Rs + 50))
  ETP_D = rnorm(nsynth, ETP_D_m, 0.2)  #PET in mm/day, Turc 1961
  ETP_M <- ETP_D * 30  #mm/month
  
  #Actual ET
  ETA_var = rnorm(nsynth, 1, 0.2) #This noise parmeter limits ETA<CMP_mm but allows variation around ETP, as observed
  ETA = CMP_mm*3 * (1 / (sqrt(1 + (1 / ((ETP_M / (CMP_mm*3)) * ETA_var)) ^ 2))) #AET in mm/month from Budyko curve
  #here scaled eta to quarter precip, assuming potential carry-over
  
  #Free air porosity
  #Have updated, now scales volumetrically w/ excess precipitation relative to pore space
  FAP <- pmin((pores - (CMP_mm - ETA)/(L*10*pores)), pores)
  FAP = pmax(FAP,0.01) #dimensionless
  
  #CO2 Diffusion coefficients
  DIFC = FAP * tort * 0.1369 * (CQT_K / 273.15) ^ 1.958
  
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
  
  #Soil carbonate C isotopes
  A_CO2_Carb <- 2.71828 ^ (-2.988e3 / CQT_K ^ 2 + 7.6663 / CQT_K - 0.0024612)
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
  DIFO <- 1.637e-8 * (CQT_K / 216.25 - 1) ^ 2.074 * (pores - FAP) * tort   ## should be soil water fraction, 
  ## pores - FAP. units: m2/sec. However, the the paper assumes total saturation, where FAP = 0
  z_i <- DIFO / E_s #mean penetration depth of evap, in m
  
  #Soil water O isotopes
  DRF <- 1 + 0.8 * (1 / 0.9723 - 1)
  R_O_surface <- ((1 - h) * DRF * R_O_P + h * R_O_atm) / (1 / A_atmP)
  R_O_soil <- ((R_O_surface - R_O_P) * 2.71828 ^ (-z_m / z_i)) + R_O_P
  R_O_soil = R_O_soil * esw + R_O_P * (1 - esw)  #soil water is esw % evaporated fraction
  dO_soil <- (R_O_soil/RO.vsmow - 1) * 1000
  
  #Soil carbonate O isotopes
  A_O <- 2.71828 ^ ((2.78e6 / CQT_K ^ 2 - 2.89) / 1000)
  R_O_Carb <- R_O_soil * A_O
  
  dC_Carb <- (R_Carb / RC.vpdb - 1) * 1000
  dO_Carb <- (R_O_Carb / RO.vpdb - 1) * 1000
  
  dat = c(median(z), median(dO_soil), median(dO_P), median(dC_Carb), median(dO_Carb))
  
  return(dat)
}

# This code was used for testing and making some validation plots...

nsynth=5000
sites = read.csv("valsites.sel.csv")
data.comp = merge.data.frame(sites, data.aves, by.x = "Site", by.y = "Group.1")

hq_pred = data.frame(d13C=numeric(0), d18O=numeric(0))
for(i in 1: nrow(sites)){
  hq_pred[i,] = sm_forward(data.comp$map.wc[i], data.comp$mat.wc[i], data.comp$hqp.frac[i], data.comp$hqt.offset[i], 280)
}
hq_pred$Site = data.comp$Site

hq.comp = merge.data.frame(hq_pred, data.comp, by.x = "Site", by.y = "Site", all.x=TRUE)
dq.comp = merge.data.frame(dq_pred, data.comp, by.x = "Site", by.y = "Site", all.x=TRUE)

c = ceiling((hq.comp$map.wc) / max(hq.comp$map.wc) * 6)
pal = rainbow(6)

jpeg("validation.jpg", res=300, units="in", width = 10, height = 5)
plot(hq.comp$d18O, hq.comp$d18O.measured, pch=16, col=pal[c],
     xlab=expression(paste("Predicted ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Observed ",delta^{18}, "O (\u2030)")))
points(hq.comp$d18O[hq.comp$map.wc>100], hq.comp$d18O.measured[hq.comp$map.wc>100], pch=1)
abline(0,1)
dev.off()

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
par(mar=c(5,5,1,1))
plot(hq.comp$d13C, hq.comp$d13C.measured, pch=16, col=pal[c], 
     xlab=expression(paste("Predicted ",delta^{13}, "C (\u2030)")),
     ylab=expression(paste("Observed ",delta^{13}, "C (\u2030)")))
points(hq.comp$d13C[hq.comp$map.wc>100], hq.comp$d13C.measured[hq.comp$map.wc>100], pch=1)
abline(0,1)
legend("bottomright", title = "MAP (mm)", fill=rainbow(6), legend=c("100 - 122", "122 - 244", "244 - 366" , "366 - 488", "488 - 600", "600 - 722"))

plot(hq.comp$d18O, hq.comp$d18O.measured, pch=16, col=pal[c],
     xlab=expression(paste("Predicted ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Observed ",delta^{18}, "O (\u2030)")))
points(hq.comp$d18O[hq.comp$map.wc>100], hq.comp$d18O.measured[hq.comp$map.wc>100], pch=1)
abline(0,1)

par(mar=c(5,5,1,1))
plot(dq.comp$d13C, dq.comp$d13C.measured, pch=16, col=pal[c], 
     xlab=expression(paste("Predicted ",delta^{13}, "C (\u2030)")),
     ylab=expression(paste("Observed ",delta^{13}, "C (\u2030)")))
points(dq.comp$d13C[dq.comp$map.wc>100], dq.comp$d13C.measured[dq.comp$map.wc>100], pch=1)
abline(0,1)

plot(dq.comp$d18O, dq.comp$d18O.measured, pch=16, col=pal[c],
     xlab=expression(paste("Predicted ",delta^{18}, "O (\u2030)")),
     ylab=expression(paste("Observed ",delta^{18}, "O (\u2030)")))
points(dq.comp$d18O[dq.comp$map.wc>100], dq.comp$d18O.measured[dq.comp$map.wc>100], pch=1)
abline(0,1)

###############Forward model function for use in sensitivity testing
esw = 1
spre = 0
sm_optimizer = function(MAP, MAT, P_seas, T_seas, pCO2, spre, esw){
  
  deltaA = rnorm(nsynth, -6.5, 0.3)
  pores = rnorm(nsynth, 0.46, 0.1)
  tort = rnorm(nsynth, 0.6, 0.1)
  
  #Solar radiation, here fixed
  Rs = 20.35
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
  
  #Convert pCO2 to units mol/cm^3
  pCO2_mcc = pCO2 / (0.08206 * CQT_K * 10^9)  #mol/cm^3
  
  
  #Relative humidity, now beta distribution
  h_m <- 0.25 + 0.7 * (CQP / 900)
  h_var = 0.05^2
  size = h_m*(1-h_m)/h_var - 1
  alpha = h_m * size
  beta = (1-h_m) * size
  h = rbeta(nsynth, alpha, beta)
  RH <- h * 100
  
  #Precipitation O isotope ratios 
  dO_P_m <- -13.7 + 0.55 * (MAT + T_seas * spre)  #Relevant precip is spre % from CQ
  dO_P = rnorm(nsynth, dO_P_m, 1.7)
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
  
  #Respiration rate, now gamma dist
  R_month_m <- 1.25 * exp(0.055 * CQT) * CMP_cm / (4.78 + CMP_cm)  #Raich 2002, gC/m2day
  theta = (R_month_m*0.5)^2/R_month_m #gamma scale parameter, using mean residual of 50% based on Raich validation data 
  k = R_month_m / theta #gamma shape parameter
  R_month = rgamma(nsynth, shape = k, scale = theta) #lets use gamma for these quants bounded at zero....
  R_month = R_month / (12.01 * 100^2)  #molC/cm2month
  R_sec <- R_month / (24 * 3600)  #molC/cm2s
  R_sec = R_sec / L / pores #molC/cm3s
  
  #Potential ET
  ETP_D_m <- ifelse (RH < 50, 0.0133 * (CQT / (CQT + 15)) * (1/23.885 * Rs + 50) * (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (1/23.885 * Rs + 50))
  ETP_D = rnorm(nsynth, ETP_D_m, 0.2)  #PET in mm/day, Turc 1961
  ETP_M <- ETP_D * 30  #mm/month
  
  #Actual ET
  ETA_var = rnorm(nsynth, 1, 0.2) #This noise parmeter limits ETA<CMP_mm but allows variation around ETP, as observed
  ETA = CMP_mm*3 * (1 / (sqrt(1 + (1 / ((ETP_M / (CMP_mm*3)) * ETA_var)) ^ 2))) #AET in mm/month from Budyko curve
  #here scaled eta to quarter precip, assuming potential carry-over
  
  #Free air porosity
  #Have updated, now scales volumetrically w/ excess precipitation relative to pore space
  FAP <- pmin((pores - (CMP_mm - ETA)/(L*10*pores)), pores)
  FAP = pmax(FAP,0.01) #dimensionless
  
  #CO2 Diffusion coefficient
  DIFC = FAP * tort * 0.1369 * (CQT_K / 273.15) ^ 1.958
  
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
  
  #Soil carbonate C isotopes
  A_CO2_Carb <- 2.71828 ^ (-2.988e3 / CQT_K ^ 2 + 7.6663 / CQT_K - 0.0024612)
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
  DIFO <- 1.637e-8 * (CQT_K / 216.25 - 1) ^ 2.074 * (pores) * tort   ## should be soil water fraction, 
  ## pores - FAP. units: m2/sec. However, the the paper assumes total saturation, where FAP = 0. 
  ## pores - FAP gives no evap at low precip regimes bc DIFO is almost 0, which makes the model insensitive to ews
  z_i <- DIFO / E_s #mean penetration depth of evap, in m
  
  #Soil water O isotopes
  DRF <- 1 + 0.8 * (1 / 0.9723 - 1)
  R_O_surface <- ((1 - h) * DRF * R_O_P + h * R_O_atm) / (1 / A_atmP)
  R_O_soil <- ((R_O_surface - R_O_P) * 2.71828 ^ (-z_m / z_i)) + R_O_P
  R_O_soil = R_O_soil * esw + R_O_P * (1 - esw)  #soil water is esw % evaporated fraction
  dO_soil <- (R_O_soil/RO.vsmow - 1) * 1000
  
  #Soil carbonate O isotopes
  A_O <- 2.71828 ^ ((2.78e6 / CQT_K ^ 2 - 2.89) / 1000)
  R_O_Carb <- R_O_soil * A_O
  
  dC_Carb <- (R_Carb / RC.vpdb - 1) * 1000
  dO_Carb <- (R_O_Carb / RO.vpdb - 1) * 1000
  
  return(median(dO_Carb))
}

sites = sites[sites$map.wc > 100,]

s = seq(0,0.95,0.05)
ss = seq(0, 1-0.05/20, 0.05/20)
ss = ss*20
ss = trunc(ss)
ss = ss/20

parms = data.frame(spres = ss, esws = rep(s, 20), rmse = numeric(400))

for(j in 1:nrow(parms)){
  opt = numeric()
  for(i in 1: nrow(sites)){
    opt[i] = sm_optimizer(sites$map.wc[i], sites$mat.wc[i], sites$hqp.frac[i], sites$hqt.offset[i], 280, parms$spres[j], parms$esws[j])
  }

  opt = data.frame(Site = sites$Site, d18O = opt)
  opt = merge.data.frame(opt, data.comp, by.x = "Site", by.y = "Site", all.x=TRUE)
  
  mse = (opt$d18O - opt$d18O.measured)^2
  mse = mean(mse)
  parms$rmse[j] = sqrt(mse)
}

View(parms)

rmses = matrix(parms$rmse, 20, 20)
rmses = rmses[c(20:1),]
rmses.rast = raster(rmses, xmn=0, xmx=0.95, ymn=0, ymx=0.95)

par(mai=c(1.05, 0.9, 0.6, 0.8))
plot(rmses.rast, xlab="% seasonal rainfall", ylab="% evaporated water")  #now need to make a nice plot...
mtext("RMSE, per mil", 4, line=1.9)
dev.off()
jpeg("O_opt.jpg", res=300, units="in", width = 5.7, height = 5)
## Repeat for DQ
mean(sites$dqt.offset)
for(j in 1:nrow(parms)){
  opt = numeric()
  for(i in 1: nrow(sites)){
    opt[i] = sm_optimizer(sites$map.wc[i], sites$mat.wc[i], sites$dqp.frac[i], sites$dqt.offset[i], 280, parms$spres[j], parms$esws[j])
  }
  
  opt = data.frame(Site = sites$Site, d18O = opt)
  opt = merge.data.frame(opt, data.comp, by.x = "Site", by.y = "Site", all.x=TRUE)
  
  mse = (opt$d18O - opt$d18O.measured)^2
  mse = mean(mse)
  parms$rmse[j] = sqrt(mse)
}

View(parms)

rmses = matrix(parms$rmse, 20, 20)
rmses = rmses[c(20:1),]
rmses.rast = raster(rmses, xmn=0, xmx=0.95, ymn=0, ymx=0.95)

jpeg("O_opt.jpg", res=300, units="in", width = 5.7, height = 5)
par(mai=c(1.05, 0.9, 0.6, 0.8))
plot(rmses.rast, xlab="% seasonal rainfall", ylab="% evaporated water")  #now need to make a nice plot...
mtext("RMSE, per mil", 4, line=1.9)
dev.off()


### Optimizing the model to fraction of estimated respiration rate (hot quarter)

sm_optimizer_r = function(MAP, MAT, P_seas, T_seas, pCO2, rr){
  
  deltaA = rnorm(nsynth, -6.5, 0.3)
  pores = rnorm(nsynth, 0.46, 0.1)
  tort = rnorm(nsynth, 0.6, 0.1)
  
  #Solar radiation, here fixed
  Rs = 20.35
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
  
  #Convert pCO2 to units mol/cm^3
  pCO2_mcc = pCO2 / (0.08206 * CQT_K * 10^9)  #mol/cm^3
  
  
  #Relative humidity, now beta distribution
  h_m <- 0.25 + 0.7 * (CQP / 900)
  h_var = 0.05^2
  size = h_m*(1-h_m)/h_var - 1
  alpha = h_m * size
  beta = (1-h_m) * size
  h = rbeta(nsynth, alpha, beta)
  RH <- h * 100
  
  #Precipitation O isotope ratios 
  dO_P_m <- -13.7 + 0.55 * (MAT + T_seas)
  dO_P = rnorm(nsynth, dO_P_m, 1.7)
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
  
  #Respiration rate, now gamma dist
  R_month_m <- rr* 1.25 * exp(0.055 * CQT) * CMP_cm / (4.78 + CMP_cm)  #Raich 2002, gC/m2day
  theta = (R_month_m*0.5)^2/R_month_m #gamma scale parameter, using mean residual of 50% based on Raich validation data 
  k = R_month_m / theta #gamma shape parameter
  R_month = rgamma(nsynth, shape = k, scale = theta) #lets use gamma for these quants bounded at zero....
  R_month = R_month / (12.01 * 100^2)  #molC/cm2month
  R_sec <- R_month / (24 * 3600)  #molC/cm2s
  R_sec = R_sec / L / pores #molC/cm3s
  
  #Potential ET
  ETP_D_m <- ifelse (RH < 50, 0.0133 * (CQT / (CQT + 15)) * (1/23.885 * Rs + 50) * (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (1/23.885 * Rs + 50))
  ETP_D = rnorm(nsynth, ETP_D_m, 0.2)  #PET in mm/day, Turc 1961
  ETP_M <- ETP_D * 30  #mm/month
  
  #Actual ET
  ETA_var = rnorm(nsynth, 1, 0.2) #This noise parmeter limits ETA<CMP_mm but allows variation around ETP, as observed
  ETA = CMP_mm*3 * (1 / (sqrt(1 + (1 / ((ETP_M / (CMP_mm*3)) * ETA_var)) ^ 2))) #AET in mm/month from Budyko curve
  #here scaled eta to quarter precip, assuming potential carry-over
  
  #Free air porosity
  #Have updated, now scales volumetrically w/ excess precipitation relative to pore space
  FAP <- pmin((pores - (CMP_mm - ETA)/(L*10*pores)), pores)
  FAP = pmax(FAP,0.01) #dimensionless
  
  #CO2 Diffusion coefficient
  DIFC = FAP * tort * 0.1369 * (CQT_K / 273.15) ^ 1.958
  
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
  
  #Soil carbonate C isotopes
  A_CO2_Carb <- 2.71828 ^ (-2.988e3 / CQT_K ^ 2 + 7.6663 / CQT_K - 0.0024612)
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
  DIFO <- 1.637e-8 * (CQT_K / 216.25 - 1) ^ 2.074 * pores * tort  ##check these, why the 2 lines?
  z_i <- DIFO / E_s #mean penitration depth of evap, in m
  
  #Soil water O isotopes
  A_O <- 2.71828 ^ ((2.78e6 / CQT_K ^ 2 - 2.89) / 1000)
  R_O_Carb <- R_O_P * A_O
  
  dC_Carb <- (R_Carb / RC.vpdb - 1) * 1000
  dO_Carb <- (R_O_Carb / RO.vpdb - 1) * 1000
  
  return(median(dC_Carb))
}

sites = sites[sites$map.wc > 100,]

rr = seq(0.01, 1, 0.01)
rr

parms = data.frame(rr = rr, rmse = numeric(100))

## All sites - Hot quarter

for(j in 1:nrow(parms)){
  opt = numeric()
  for(i in 1: nrow(sites)){
    opt[i] = sm_optimizer_r(sites$map.wc[i], sites$mat.wc[i], sites$hqp.frac[i], sites$hqt.offset[i], 280, parms$rr[j])
  }
  
  opt = data.frame(Site = sites$Site, d13C = opt)
  opt = merge.data.frame(opt, data.comp, by.x = "Site", by.y = "Site", all.x=TRUE)
  
  mse = (opt$d13C - opt$d13C.measured)^2
  mse = mean(mse)
  parms$rmse[j] = sqrt(mse)
}

View(parms)
plot(parms$rmse ~ parms$rr, type="n",main = "Respiration Optimization Hot Quarter", ylab = "RMSE", xlab = "Fraction of Estimated Respiration")
lines(parms$rmse ~ parms$rr)
## Dry Quarter

for(j in 1:nrow(parms)){
  opt = numeric()
  for(i in 1: nrow(sites)){
    opt[i] = sm_optimizer_r(sites$map.wc[i], sites$mat.wc[i], sites$dqp.frac[i], sites$dqt.offset[i], 280, parms$rr[j])
  }
  
  opt = data.frame(Site = sites$Site, d13C = opt)
  opt = merge.data.frame(opt, data.comp, by.x = "Site", by.y = "Site", all.x=TRUE)
  
  mse = (opt$d13C - opt$d13C.measured)^2
  mse = mean(mse)
  parms$rmse[j] = sqrt(mse)
}

View(parms)
plot(parms$rmse ~ parms$rr, type="n", main = "Respiration Optimization Dry Quarter", ylab = "RMSE", xlab = "Fraction of Estimated Respiration")
lines(parms$rmse ~ parms$rr)
## Only clumped sites

parms = data.frame(rr = rr, rmse = numeric(400))

for(j in 1:nrow(parms)){
  opt = numeric()
  for(i in 1: nrow(clump.comp.meantemp)){
    opt[i] = sm_optimizer_r(clump.comp.meantemp$map.wc[i], clump.comp.meantemp$mat.wc[i], clump.comp.meantemp$Clumped.frac[i], clump.comp.meantemp$Informed.offset[i], 280, parms$rr[j])
  }
  
  opt = data.frame(Site = clump.comp.meantemp$Site, d13C = opt)
  opt = merge.data.frame(opt, data.comp, by.x = "Site", by.y = "Site", all.x=TRUE)
  
  mse = (opt$d13C - opt$d13C.measured)^2
  mse = mean(mse)
  parms$rmse[j] = sqrt(mse)
}

View(parms)
plot(parms$rr ~ parms$rmse, main = "Clumped Data Respiration Optimization", ylab = "Fraction of Estimated Respiration", xlab = "RMSE")
