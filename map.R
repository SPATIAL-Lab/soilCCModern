## Make a map of sites used
library(maps)
library(scales)

# Lat - long of sites is already in site data

jpeg("map.jpg", res=300, width=10, height=5, units="in")
map("world", fill=T, col=grey(0.6, alpha = 0.7), bg="white", mar=c(0,0,0,0), xlim=c(-180, 180), ylim=c(-90,90), lty=0)
points(sites$Lon, sites$Lat, col=alpha("red", alpha=0.6), pch=16, cex=1, lwd = 1)
dev.off()

