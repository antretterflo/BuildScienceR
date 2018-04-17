# Hygrothermal Functions --------------------------------------------------

# Functions to compute hygrothermal values
# Florian Antretter
# Date: 2018-04-12
# Version: 0.1


# Functions ---------------------------------------------------------------

# calcFtoC : Compute temperature in Celsius from Fahrenheit
calcFtoC <- function(Temperature_F) {
  (Temperature_F - 32) * (5/9)
}

# calcCtoF : Compute temperature in Fahrenheit from Celsius
calcCtoF <- function(Temperature_C) {
  Temperature_C * (9/5) + 32
}


# Conversion formulas according to Buck -----------------------------------

# Original reference is
# Buck Research Instruments, LLC (May 2012):
# Model CR-1A Hygrometer with autofill - operating manual
# Appendix 1: Humidity conversion equations (according to A. Buck 1981)

# Td: Dew Point Temperature in C
# Tc: Temperature in C
# Tk: Absolute Temperature in K
# AP: Absolute Pressure in mbar
# VP: Vapor Pressure in mbar
# VPs: Saturation Vapor Pressure in mbar
# MR: Mixing Ratio by weight in ppm
# AH: Absolute Humidity in g/m3
# AHs: Saturation Absolute Humidity in g/m3
# RH: Relative Humidity in %
# Ent: Enthalpy in J/kg

### Function enhancement factor
EF <- function(Tc, AP = 1013) {
  ifelse(Tc < 0, 
         (1 + 10^-4 * (2.2 + AP * (0.0383 + 6.4 * 10^-6 * Tc^2))), 
         (1 + 10^-4 * (7.2 + AP * (0.0320 + 5.9 * 10^-6 * Tc^2))))
}

### Function 1: VPs vs. T
f1 <- function(Tc, AP = 1013) {
  con1 <- ifelse(Tc < 0, 6.1115, 6.1121)
  con2 <- ifelse(Tc < 0, 23.036, 18.678)
  con3 <- ifelse(Tc < 0, 333.7, 234.5)
  con4 <- ifelse(Tc < 0, 279.82, 257.14)
  EF(Tc, AP) * con1 * exp((con2 - (Tc / con3)) * (Tc / (con4 + Tc)))
}

### Function 2: Td vs. VP
f2 <- function(Tc, RH, AP = 1013) {
  con1 <- ifelse(Tc < 0, 6.1115, 6.1121)
  con2 <- ifelse(Tc < 0, 23.036, 18.678)
  con3 <- ifelse(Tc < 0, 333.7, 234.5)
  con4 <- ifelse(Tc < 0, 279.82, 257.14)
  
  s <- log(((RH / 100) * f1(Tc, AP)) / EF(Tc, AP)) - log(con1)
  
  (con3 / 2) * (con2 - s - ((con2 - s)^2 - 4 * con4 * (s / con3))^0.5)
}

### Temperature dependent values
# Vapor Pressure in Pa
calcVPs_Pa <- function(Tc, AP = 1013) {
  f1(Tc, AP) * 100
}

# Absolute Humidity in kg/m3
calcAHs_kgm3 <- function(Tc, AP = 1013) {
  (0.2167 * f1(Tc, AP) / (Tc + 273.15))
}

# Mixing Ratio in kg/kg
calcMRs_kgkg <- function(Tc, AP = 1013) {
  VP <- f1(Tc, AP)
  (0.622 * VP) / (AP - VP)
}

# Dew Point Temperature in C
calcTd_C <- function(Tc, RH, AP = 1013) {
  f2(Tc, RH, AP)
}

# Enthalpy in J/kg 
# see also https://www.engineeringtoolbox.com/enthalpy-moist-air-d_683.html
# or https://www.wikihow.com/Calculate-the-Enthalpy-of-Moist-Air
calcEnt_Jkg <- function(Tc, RH, AP = 1013) {
  specHeatAir <- 1006
  specHeatWV0 <- 2501000
  specHeatWV  <- 1840
  specHeatAir * Tc + ((RH / 100) * calcMRs_kgkg(Tc, AP)) * (specHeatWV0 + specHeatWV * Tc)
}
