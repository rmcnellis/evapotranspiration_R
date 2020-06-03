# evapotranspiration_R
This repository contains the R functions necessary for calculating evapotranspiration in R as in McNellis & Smith Winter cover crops increase winter albedo and latent heat flux in a Texas High Plains agro-ecosystem.

The script calc_ET.R contains code for two evapotranspiration functions. Both functions are dependent on the packages plantecophys and bigleaf in order to calculate wet bulb temperature and vapor pressure deficit. You must run 'install.packages(“plantecophys”)' and 'install.packages(“bigleaf”)' prior to running the code.

## DOI badge
[![DOI](https://zenodo.org/badge/264007622.svg)](https://zenodo.org/badge/latestdoi/264007622)

## calc_dual_ET
This function uses the dual-source method for calculating evapotranspiration adapted from Shuttleworth & Wallace 1985. The function will calculate:
- net radiation
- ground heat flux 
- bulk stomatal resistance of the canopy
- bulk boundary layer resistance of canopy
- aerodynamic resistance between the substrate and canopy
- aerodynamic resistance between canopy and boundary layer
- canopy and surface coefficients, canopy and surface evapotranspiration
- canopy plus surface evapotranspiration

The inputs required are:
- maximum, minimum, and average air temperature
- air pressure
- relative humidity
- day of year
- latitude
- elevation
- sunshine duration in a day
- wind speed
- albedo
- leaf area index
- canopy height
- stomatal conductance
- soil moisture
- saturation soil moisture
- soil surface roughness
- psychrometer constant
- fraction of extraterrestrial radiation reaching the earth on a clear day
- Schmidt and Prandtl numbers
- reference height for temperature and wind speed
- eddy diffusivity decay constant in a crop with complete canopy cover
- specific heat of air

## calc_soil_ET
This function uses the Penman-Monteith equation for calculating soil evapotranspiration with no canopy cover. The function will calculate:
- net radiation
- ground heat flux
- soil surface resistance
- aerodynamic resistance
- evapotranspiration

The inputs required are:
- maximum, minimum, and average air temperature
- air pressure
- relative humidity
- day of year
- latitude
- elevation
- sunshine duration in a day
- wind speed
- albedo
- soil moisture
- saturation soil moisture
- soil surface roughness
- psychrometer constant
- fraction of extraterrestrial radiation reaching the earth on a clear day
- Schmidt and Prandtl numbers
- reference height for temperature and wind speed
- eddy diffusivity decay constant in a crop with complete canopy cover
- specific heat of air

The script test_calc_ET.R will test the functions. It is suggested that users run this script first to ensure the functions will run properly.

Any questions or issues can be submitted via GitHub or directed to Risa McNellis (risa.mcnellis@gmail.com).
