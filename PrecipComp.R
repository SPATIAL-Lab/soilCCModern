precipcomp <- merge.data.frame(precipcomp, sites, by.x="Site", by.y="Site")

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