## Make a map of sites ued

# Lat - long of sites is already in site data

sites$Lat
sites$Lon
jpeg("map.jpg", res=300, width=10, height=5, units="in")
map("world", fill=T, col="white", bg="lightblue", mar=c(0,0,0,0))
points(sites$Lon, sites$Lat, col=alpha("red", 0.65), pch=16, cex=1.5)
dev.off()