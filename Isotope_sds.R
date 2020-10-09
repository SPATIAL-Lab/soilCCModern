# Site isotope sds

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