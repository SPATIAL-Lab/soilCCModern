## Respiration data subset into subhumid to arid

srdb = read.csv("srdb-data.csv")
srdb_carb <- subset(srdb, srdb$MAP < 300)
