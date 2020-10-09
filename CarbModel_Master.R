###Code for conducting validation of soil carbonate model
#load required packages
library(scales)
library(readxl)
library(RColorBrewer)

# Read in input data R workspace
load("Input_Data_Final.RData")

#Inter-site typical isotope standard deviations
source("Isotope_sds.R", echo=T)

#Code calculating d18O-temperature relationship
source("T_O_corr.R", echo=T)

#Run theoretical model
source("Theoretical.R", echo=T)

#Run theoretical model with depth as input
source("Theoretical_depth.R", echo=T)

#Run theoretical model with OIPC precipitation
source("Theoretical_OIPC.R", echo=T)

#Run Seasonal Precip and Evap Optimization
source("Evap_Seasonal_OPT.R", echo=T)

#Run Seasonal Precip and Evap Optimization with depth as input
source("Evap_Seasonal_OPT_depth.R", echo=T)

#Run respiration optimization
source("Respiration_OPT.R", echo=T)

#Run respiration optmiziation with depth as input
source("Respiration_OPT_depth.R", echo=T)

#Apply optimization
source("Optimized.R", echo=T)

#Apply optimization with depth as input
source("Optimized_depth.R", echo=T)

#Sensitivity tests
source("SensTest.R", echo=T)

#Prediction maps - worldwide
source("PredMaps_OPT.R", echo=T)
