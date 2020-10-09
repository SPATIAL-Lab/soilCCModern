## README for soilCCModern, code repository for Fischer-Femal and Bowen (2020) 

# Brief description of R scripts and data inputs

# Optimized_Carb_Model.R has the optimized model for both warm and dry quarter conditions without anything else. This script can be used to predict d18O and d13C carbonate for any climate conditions.

# CarbModel_Master.R is the master R code for loading the data and running all other code
	# Isotope_sds.R calculates typical standard deviations of carbonate isotope values at a site
	# T_O_corr.R runs statistics on the d18Op-T correlation and OIPC d18Op values
	# Theoretical.R runs the theoretical model predictions and comparison to modern data
	# Theoretical_depth.R runs the measured-depth theoretical model predictions and comparison to modern data
	# Theoretical_OIPC.R runs the theoretical model with OIPC oxygen isotope values instead of using the d18O-T relationship
	# Respiration_OPT.R runs the optimization of respiration rate 
	# Respiration_OPT_depth.R runs the optimization of respiration rate for the measured-depth model
	# Evap_Seasonal_OPT.R runs the optimization for seasonal precipitation and evaporative effects 
	# Evap_Seasonal_OPT_depth.R runs the optimization for seasonal precipitation and evaporative effects for the measured-depth model
	# Optimized.R runs the optimized model and comparison to modern data
	# Optimized_depth.R runs the optimized measured-depth model and comparison to modern data
	# SensTest.R runs the sensitivity tests
	# PredMaps_OPT.R runs the spatial predictions of the optimized model

# Input_Data_Final.RData is the input carbonate data and the WorldClim climate variables at each site
# rasters.RData contains the WorldClim rasters needed to use PredMaps_OPT.R and predict carbonate isotope values worldwide