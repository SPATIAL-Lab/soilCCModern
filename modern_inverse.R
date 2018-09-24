library(R2OpenBUGS)
library(coda)
library(rjags)
library(R2jags)

set.seed(1234)

soil_model = " model {

  dO_Carb ~ dnorm(dO_Carb_m, 1 / stdevO ^ 2)
  dC_Carb ~ dnorm(dC_Carb_m, 1 / stdevC ^ 2)

  dO_Carb_m <- (R_O_Carb / RO.vpdb - 1) * 1000
  dC_Carb_m <- (R_Carb / RC.vpdb - 1) * 1000

  R_O_Carb <- R_O_P * A_O
  A_O <- 2.71828 ^ ((2.78e6 / CQT_K ^ 2 - 2.89) / 1000)

  R_Carb <- R_Soil / A_CO2_Carb
  R_Soil <- (dC_Soil / 1000 + 1) * RC.vpdb
  A_CO2_Carb <- 2.71828 ^ (-2.988e3 / CQT_K ^ 2 + 7.6663 / CQT_K - 0.0024612)

  dC_Soil <- (dC_Soil.num / (dC_Soil.denom * RC.vpdb) - 1) * 1000
  dC_Soil.denom <- dC_Soil.resp * (1 - DIFC * deltaP_hat / DIFC13) + pCO2_mcc * (1 - deltaA_hat)
  dC_Soil.num <- dC_Soil.resp * DIFC * deltaP_hat / DIFC13 + pCO2_mcc * deltaA_hat
  dC_Soil.resp <- R_sec/(DIFC) * (L * z - z^2 / 2)
  deltaA_hat <- (deltaA / 1000 + 1) * RC.vpdb / (1 + RC.vpdb * (deltaA / 1000 + 1))
  deltaP_hat <- (deltaP / 1000 + 1) * RC.vpdb / (1 + RC.vpdb * (deltaP / 1000 + 1))

  deltaP <- deltaA - (deltaP_pCO2 - W)

  deltaP_pCO2 ~ dnorm(deltaP_pCO2_m, 1 / 0.5 ^ 2)
  deltaP_pCO2_m <- 28.26 * 0.35 * (pCO2 + 15) / (28.26 + 0.35 * (pCO2 + 15))

  W ~ dnorm(W_m, 1 / 0.5 ^ 2)
  W_m <- 22.65 - (1.2 * (MAP + 975)) / (27.2 + 0.04 * (MAP + 975))

  DIFC13 <- DIFC * (1 / 1.004443)
  DIFC <- EPS * 0.1369 * (CQT_K / 273.15) ^ 1.958

  EPS <- ifelse(EPS.1 > 0.01, EPS.1, 0.01)
  EPS.1 <- ifelse(EPS.2 < EPSmax, EPS.2, EPSmax)
  EPS.2 <- EPSmax - (CMP_mm - ETA)/(1000*EPSmax)

  ETA <- CMP_mm*3 * (1 / (sqrt(1 + (1 / ((ETP_M / (CMP_mm*3)) * ETA_var)) ^ 2)))
  ETA_var ~ dnorm(1, 1 / 0.2^2)

  ETP_M <- ETP_D * 30 
  ETP_D ~ dnorm(ETP_D_m, 1 / 0.2^2)
  ETP_D_m <- ifelse (RH < 50, 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50) * (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50))

  R_sec <- R_month / (24 * 3600 * L)
  R_month ~ dgamma(k.r, 1 / theta.r)
  k.r <- R_month_m / theta.r
  theta.r <- (R_month_m * 0.5)^2 / R_month_m 
  R_month_m <- 1.25 * exp(0.05452 * CQT) * CMP_cm / ((12.7 + CMP_cm) * 12.01 * 100^2)

  z ~ dgamma(k.z, 1 / theta.z)T(,100)
  k.z <- z_mean.thick / theta.z
  theta.z <- 20^2 / z_mean.thick
  z_mean.thick <- z_mean + z_thick/2
  z_thick <- abs(CMP_mm - MAP/12) * 0.74 + 17.3 
  z_mean <- MAP * 0.0925 + 13.4

  R_O_P <- (dO_P / 1000 + 1) * RO.vsmow
  dO_P ~ dnorm(dO_P_m, 1 / 1.7^2)
  dO_P_m <- -13.7 + 0.55 * MAT 

  RH <- h * 100
  h ~ dbeta(alpha.h, beta.h)
  alpha.h <- h_m * size.h
  beta.h <- (1-h_m) * size.h
  size.h <- h_m*(1-h_m)/h_var - 1
  h_var <- 0.05^2
  h_m <- 0.25 + 0.7 * (CQP / 900)

  pCO2_mcc <- pCO2 / (v.million * avo)
  v.million <- 0.0821e9 * CQT_K / avo  
  avo <- 6.023e23

  CQT_K <- CQT + 273
  CQT <- MAT + T_seas
  CMP_cm <- CQP / 30
  CMP_mm <- CQP / 3
  CQP <- ifelse(MAP * P_seas < 3, 3, MAP * P_seas) 

  MAP ~ dunif(50, 1000)
  P_seas ~ dbeta(3, 12)
  MAT ~ dnorm(12, 1 / 10 ^ 2)
  T_seas ~ dnorm(10, 1 / 5 ^ 2)

  EPSmax <- 0.1 + ST * 0.4
  ST ~ dnorm(0.5, 1 / 0.05 ^ 2)I(0, 1)

  deltaA ~ dnorm(-6.5, 1 / 0.1 ^ 2)
  pCO2 ~ dnorm(280, 1 / 10^2)

  Rs <- 20.35
  RC.vpdb = 0.011237
  RO.vsmow = 0.0020052
  RO.vpdb = 0.002067181
  L = 100

} "

smold = " model {

dO_Carb ~ dnorm(dO_Carb_m, 1 / stdevO ^ 2)
dC_Carb ~ dnorm(dC_Carb_m, 1 / stdevC ^ 2)

dO_Carb_m <- (R_O_Carb / 0.0020672 - 1) * 1000
dC_Carb_m <- (R_Carb / 0.011237 - 1) * 1000


  R_O_Carb <- R_O_P * A_O
  A_O <- 2.71828 ^ ((2.78e6 / CQT_K ^ 2 - 2.89) / 1000)

  R_Carb <- R_Soil / A_CO2_Carb
  R_Soil <- (dC_Soil / 1000 + 1) * RC.vpdb
  A_CO2_Carb <- 2.71828 ^ (-2.988e3 / CQT_K ^ 2 + 7.6663 / CQT_K - 0.0024612)

  dC_Soil <- (dC_Soil.num / (dC_Soil.denom * RC.vpdb) - 1) * 1000
  dC_Soil.denom <- dC_Soil.resp * (1 - DIFC * deltaP_hat / DIFC13) + pCO2_mcc * (1 - deltaA_hat)
  dC_Soil.num <- dC_Soil.resp * DIFC * deltaP_hat / DIFC13 + pCO2_mcc * deltaA_hat
  dC_Soil.resp <- R_sec/(DIFC) * (L * z - z^2 / 2)
  deltaA_hat <- (deltaA / 1000 + 1) * RC.vpdb / (1 + RC.vpdb * (deltaA / 1000 + 1))
  deltaP_hat <- (deltaP / 1000 + 1) * RC.vpdb / (1 + RC.vpdb * (deltaP / 1000 + 1))
  R_sec <- R_month / (24 * 3600 * L)

deltaP <- deltaA - (deltaP_pCO2 - W)
deltaP_pCO2 ~ dnorm(deltaP_pCO2_m, 1 / 0.5 ^ 2)
deltaP_pCO2_m <- 28.26 * 0.35 * (pCO2 + 15) / (28.26 + 0.35 * (pCO2 + 15))

W ~ dnorm(W_m, 1 / 0.5 ^ 2)
W_m <- 22.65 - (1.2 * (MAP + 975)) / (27.2 + 0.04 * (MAP + 975))

DIFC13 <- DIFC * (1 / 1.004443)
DIFC <- EPS * p * 0.14
EPS <- ifelse(EPSmax * (ETA / CMP_mm) > EPSmax, EPSmax, EPSmax * (ETA / CMP_mm))

ETA ~ dnorm(ETA_m, 1 / 3 ^ 2)
ETA_m <- CMP_mm * (1 / (sqrt(1 + (1 / (ETP_M / CMP_mm)) ^ 2)))

ETP_M <- ETP_D * 30

ETP_D ~ dnorm(ETP_D_m, 1 / 0.2 ^ 2)
ETP_D_m <- ifelse (RH < 50, 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50) * (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50))

A <- 1/(2.71828^((-2.988e6 / CQT_K ^ 2 + 7.6663e3 / CQT_K - 2.4612) / 1000))

  R_month ~ dgamma(k.r, 1 / theta.r)
  k.r <- R_month_m / theta.r
theta.r <- (R_month_m * 0.5)^2 / R_month_m 
R_month_m <- 1.25 * exp(0.05452 * CQT) * CMP_cm / ((12.7 + CMP_cm) * 12.01 * 100^2)


z_m <- z_0 * 0.01
z_100 <- ifelse (z_0 < 100, z_0, 100) 
z_0 <- ifelse (z > 5 , z , 5)
z ~ dnorm(z_mean, 1 / 10 ^ 2)
z_mean <- (MAP - 257) / 2.775

dO_atm <- (R_O_atm / 0.0020052 - 1) * 1000
R_O_atm <- R_O_P / A_atmP
R_O_P <- (dO_P / 1000 + 1) * 0.0020052
A_atmP <- 2.71828 ^ ((5.9702e6 / CQT_K ^ 2 - 3.2801e4 / CQT_K + 52.227) / 1000)

RH <- h * 100
h <- ifelse(h_m2 < 0.95, h_m2, 0.95)
h_m2 ~ dnorm(h_m1, 1 / 0.05 ^ 2)
h_m1 <- 0.25 + 0.7 * (CQP / 900)

dO_P ~ dnorm(dO_P_m, 1)
dO_P_m <- -14.5 + 0.55 * MAT 

  pCO2_mcc <- pCO2 / (v.million * avo)
  v.million <- 0.0821e3 * CQT_K / avo  
  avo <- 6.023e23

CQT_K <- CQT + 273
CQT <- MAT + T_seas
CMP_cm <- CQP / 30
CMP_mm <- CQP / 3
CQP <- MAP * P_seas
p <- 0.5 + ST * 0.3

  EPSmax <- 0.1 + ST * 0.4
  ST ~ dnorm(0.5, 1 / 0.05 ^ 2)I(0, 1)

  MAP ~ dunif(50, 1000)
  P_seas ~ dbeta(3, 12)
  MAT ~ dnorm(12, 1 / 10 ^ 2)
  T_seas ~ dnorm(10, 1 / 5 ^ 2)

  deltaA ~ dnorm(-6.5, 1 / 0.1 ^ 2)
  pCO2 ~ dnorm(280, 1 / 10^2)

  Rs <- 20.35
  RC.vpdb = 0.011237
  RO.vsmow = 0.0020052
  RO.vpdb = 0.002067181
  L = 100

} "

####analysis

parameters <- c("pCO2", "deltaA", "MAT", "T_seas", "MAP", "P_seas", "R_month", "ST", "z")

### Pretty good case
carbonate_data = list(dC_Carb=-7.5, dO_Carb=-5.0, stdevC = 0.5, stdevO = 0.5)

# Model Fitting

bayes_fit <- jags(model.file = textConnection(soil_model), parameters.to.save = parameters, 
                  data = carbonate_data, inits = NULL, 
                  n.chains=3, n.iter = 40000, n.burnin = 10000, n.thin = 25)

post = as.data.frame(bayes_fit$BUGSoutput$sims.list)
