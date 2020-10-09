# README for soilCCModern, code repository for Fischer-Femal and Bowen (2020) 

Brief description of R scripts and data inputs

## Optimized_Carb_Model.R
Optimized model for both warm and dry quarter conditions without anything else. This script can be used to predict d18O and d13C carbonate for any climate conditions.

## CarbModel_Master.R
Master R code for loading the data and running all other code

## Isotope_sds.R
Calculates typical standard deviations of carbonate isotope values at a site

## T_O_corr.R
Runs statistics on the d18Op-T correlation and OIPC d18Op values

## Theoretical.R
Runs the theoretical model predictions and comparison to modern data
	
## Theoretical_depth.R
Runs the measured-depth theoretical model predictions and comparison to modern data
	
## Theoretical_OIPC.R 
Runs the theoretical model with OIPC oxygen isotope values instead of using the d18O-T relationship
	
## Respiration_OPT.R
Runs the optimization of respiration rate 
	
## Respiration_OPT_depth.R
Runs the optimization of respiration rate for the measured-depth model
	
## Evap_Seasonal_OPT.R
Runs the optimization for seasonal precipitation and evaporative effects 
	
## Evap_Seasonal_OPT_depth.R
Runs the optimization for seasonal precipitation and evaporative effects for the measured-depth model
	
## Optimized.R
Runs the optimized model and comparison to modern data
	
## Optimized_depth.R
Runs the optimized measured-depth model and comparison to modern data
	
## SensTest.R
Runs the sensitivity tests
	
## PredMaps_OPT.R
Runs the spatial predictions of the optimized model

## Input_Data_Final.RData
Input carbonate data and the WorldClim climate variables at each site

## rasters.RData
WorldClim rasters needed to use PredMaps_OPT.R and predict carbonate isotope values worldwide - this file is available in the Zenodo archive for this repository
