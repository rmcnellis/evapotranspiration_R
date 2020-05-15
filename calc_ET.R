# calculate dual source evapotranspiration as in McNellis & Smith, Winter cover 
# crops increase winter albedo and latent heat flux in a Texas High Plains 
# agro-ecosystem. 

###################################################################################
# variable key
###################################################################################
# alb: albedo
# as: fraction of extraterrestrial radiation reaching the earth on an overcast day
# A: total energy flux leaving the complete crop (W/m2)
# As: total energy flux leaving the substrate (W/m2)
# bs + as: fraction of extraterrestrial radiation reaching the earth on a clear day
# Cc: canopy coefficient
# Cp: specific heat of air (J kg-1 °C-1)
# Cs: surface coefficient
# dr: inverse relative distance from the earth to the sun
# D: vapor pressure deficit (kPa)
# elev: elevation above sea level (m)
# esTw: saturation vapor pressure at wet bulb temperature (kPa)
# ET: canopy plus surface evapotranspiration (W m-2)
# gl: stomatal conductance (mol m-2 s-1)
# G: soil heat flux (W/m2), arbitrarily set to 20% of Rsurf
# Gsc: solar constant, 0.082 MJ m-2 day-1
# h: canopy height (m)
# J: day of year
# k: von Karman's constant, 0.41
# l: latitude in radians
# L: leaf area index
# n: eddy diffusivity decay constant in a crop with complete canopy cover, 2.5
# ns: actual duration of sunshine in a day (hours)
# N: maximum possible sunshine duration in a day (hours)
# p: air density (kg m-3)
# pd: pressure of dry air (Pa)
# ps: saturation vapor pressure (Pa)
# pv: vapor pressure (Pa)
# Pa: air pressure (kPa)
# PMc: canopy ET (W m-2)
# PMs: surface ET (W m-2)
# Pr: Prandtl number, default to 1 from Reynold's analogy
# ra: aerodynamic resistance between canopy and boundary layer (s m-1)
# ra_i: aerodynamic resistance between canopy and boundary layer at no canopy cover
#       (s m-1)
# ra_o: aerodynamic resistance between canopy and boundary layer at full canopy 
#       cover (s m-1)
# rA: aerodynamic resistance coefficient
# rb: boundary layer resistance (s m-1)
# rc: bulk stomatal resistance of the canopy (s m-1)
# rC: canopy resistance coefficient
# rg: aerodynamic resistance between the substrate and canopy (s m-1)
# rg_i: aerodynamic resistance between the substrate and canopy at no canopy cover 
#       (s m-1)
# rg_o: aerodynamic resistance between the substrate and canopy at full canopy 
#       cover (s m-1)
# rl: leaf stomatal resistance (mol m-2 s-1)
# rs: surface resistance of the substrate (s m-1)
# rS: surface resistance coefficient
# rv: bulk boundary layer resistance of canopy (s m-1)
# Ra: extraterrestrial radiation (MJ m-2 day-1)
# Rd: specific gas constant for dry air, 287.058 J kg-1 K-1
# RH: relative humidity (%)
# Rn: net radiation (W m-2)
# Rnl: net longwave radiation (MJ m-2 day-1)
# Rns: net shortwave radiation (MJ m-2 day-1)
# Rs: shortwave radiation (MJ m-2 day-1)
# Rso: clear sky shortwave radiation (MJ m-2 day-1)
# Rsurf: radiation reaching substrate surface (W m-2)
# Rv: specific gas constant for water vapor, 461.495 J kg-1 K-1
# s: slope of saturation vapor pressure-temperature curve (kPa C-1)
# sd: solar decimation (rad)
# Sc: Schmidt number, default to 1 from Reynold's analogy
# SM: soil moisture (m3 m-3)
# SMs: saturation soil moisture (m3 m-3)
# Ta: temperature (C)
# Tmax: maximum air temperature (C)
# Tmin: minimum air temperature (C)
# Tw: wet bulb temperature (C)
# ux: wind speed at reference height (m s-1)
# Vm: volume of one mole (m3 mol-1)
# ws: sunset hour angle (rad)
# x: reference height (m)
# y: psychrometer constant (kPa C-1)
# zc: canopy roughness length of momentum (m)
# zs: soil surface roughness length of momentum (m)

# libraries
# install.packages('plantecophys')
# install.packages('bigleaf')
library(plantecophys)
library(bigleaf)


# All equation numbers refer to equations presented in McNellis & Smith

calc_dual_ET <- function(y = 0.06, Pa = 90, Ta = 25, RH = 40, J = 180, l = 0.58, elev = 900,
                         ns = 6, as = 0.25, bs = 0.5, alb = 0.2, Tmax = 30, Tmin = 20, L = 2, 
                         h = 0.3, ux = 2, x = 2, gl = 400, Sc = 1, Pr = 1, SM = 0.1, 
                         SMs = 0.5, n = 2, zs = 0.004, Cp = 1005) {
  
  # constants
  Gsc <- 0.082
  k <- 0.41
  Rd <- 287.058
  Rv <- 461.495

  # calculate environmental terms
  D <- RHtoVPD(RH = RH, TdegC = Ta, Pa = Pa)
  Tw <- wetbulb.temp(Tair = Ta, pressure = Pa * 1000, VPD = D)
  esTw <- 0.6108 * exp((17.27 * Tw) / (Tw + 238.3))
  ea <- esTw - y * Pa * (Ta - Tw)
  Vm <- (1 * 8.31441 * (Ta + 273)) / (Pa * 1000)
  s <- (4098 * (0.6108 * exp((17.27 * Ta) / (Ta + 237.3)))) / (Ta + 237.3)^2 # Eq. 18
  ps <- 6.1078 * 10^((7.5 * Ta) / (Ta + 237.3))
  pv <- ps * RH
  pd <- (Pa * 1000) - pv
  p <- (pd / (Rd * (Ta + 273.15))) + (pv / (Rv * (Ta + 273.15)))
  
  # calculate distance from the earth to the sun, solar decimation, and sunset solar angle
  dr <- 1 + 0.033 * cos(J * ((2 * pi) / 365))
  sd <- 0.409 * sin(J * ((2 * pi) / 365) - 1.39)
  ws <- (pi/2) - atan((-tan(l) * tan(sd)) / (1 - ((tan(l))^2) * ((tan(sd))^2)^0.5))
  # calculate extraterrestrial radiation
  Ra <- (24 * 60 / pi) * Gsc * dr * 
        (ws * sin(l) * sin(sd) + cos(l) * cos(sd) * sin(ws))
  N <- ws * (24 / pi)
  # calculate shortwave radiation
  Rs <- ((as + bs * (ns / N)) * Ra)
  Rso <- (0.75 + (elev * 2e-5)) * Ra
  # calculate net shortwave radiation
  Rns <- (1 - alb) * Rs
  # calculate net longwave radiation
  Rnl <- (((Tmax + Tmin) / 2) * (0.34 - 0.14 * sqrt(ea)) * (1.35 * (Rs/Rso) - 0.35))
  # calculate net radiation
  Rn_MJ <- (Rns + Rnl)
  # convert from MJ m-2 day-1 to W m-2
  Rn <- 11.82033 * Rn_MJ
  
  # calculate Rsurf, G, A, As, Eq. 22-23
  Rsurf <- Rn * exp(-0.7 * L)
  G <- (Rsurf * 0.2)
  A <- 400 - G
  As <- Rsurf - G
  
  # calculate canopy terms d, zc, uf
  d <- h * 0.63
  zc <- h * 0.13
  uf <- (0.41 * ux) / (log((x - d) / zc)) #Eq. 5
  
  # calculate resistance terms
  rl <- 1 / (gl * Vm) #Eq. 1
  rc <- rl / (2 * L) #Eq. 2
  rb <- (2 * (0.41 * uf)^-1) * ((Sc^0.67) / Pr) #Eq. 3
  rv <- rb / (2 * L) #Eq. 4
  rs <- 1 /(exp(4.2555 * (SM / SMs) - 8.206)) #Eq. 6
    # rg at no canopy cover
    rg_o <- ((log((x - d) / zc)) / ((k^2) * ux)) * 
            (h / (n * (h - d))) * (exp(n) - exp(n * (1 - ((d + zc) / h)))) # Eq. 7
    # rg at full canopy cover
    rg_i <- (log(x / zs) * log((d + zc) / zs)) / ((k^2) * ux) #Eq. 8
  # rg at given LAI
  rg <- ifelse (L <= 4, 0.25 * L * rg_o + 0.25 * (4 - L) * rg_i, rg_o) #Eq. 11-12
    # ra at no canopy cover
    ra_o <- ((log((x - d) / zc)) / ((k^2) * ux)) * 
            ((log((x - d) / (h - d))) + (h / (n * (h - d))) * 
            (exp(n * (1 - ((d + zc) / h))) - 1)) #Eq. 9
    # ra at full canopy cover
    ra_i <- ((log(x / zs))^2) / ((k^2) * ux) - rg_i #Eq. 10
  # ra at given LAI
  ra <- ifelse (L <= 4, 0.25 * L * ra_o + 0.25 * (4 - L) * ra_i, ra_o) #Eq. 11-12
  
  # calculate aerodynamic, substrate, and canopy resistance coefficients, Eq. 13-15
  rA <- (s + y) * ra 
  rS <- (s + y) * rg + (y * rs)
  rC <- (s + y) * rv + (y * rc)
  
  # calculate canopy and substrate coefficients, Eq. 16-17
  Cc <- (1 + ((rC * rA) / (rS * (rC + rA)))) ^ -1
  Cs <- (1 + ((rS * rA) / (rC * (rS + rA)))) ^ -1

  # calculate evapotranspiration using Penman-Monteith equations, Eq. 19-21
  PMc <- (s * A + (p * Cp * D - s * rv * As) / (ra + rv)) / 
         (s + y * (1 + (rc / (ra + rv))))
  PMs <- (s * A + (p * Cp * D - s * rg * (A - As)) / (ra + rg)) / 
         (s + y * (1 + rs / (ra + rg)))
  # calculate total evapotranspiration from Penman-Monteith equations and coefficients
  ET <- Cc * PMc + Cs * PMs
  
  # output
  return(data.frame(Rn = Rn, G = G, rc = rc, rv = rv, rg = rg, ra = ra, Cc = Cc, 
                    Cs = Cs, PMc = PMc, PMs = PMs, ET = ET))
}

# for bare soil only, Eq. 24
calc_soil_ET <- function(y = 0.06, Pa = 90, Ta = 25, RH = 40, J = 180, l = 0.58, elev = 900, 
                         ns = 6, as = 0.25, bs = 0.5, alb = 0.2, Tmax = 30, Tmin = 20, ux = 2, 
                         x = 2, Sc = 1, Pr = 1, SM = 0.1, SMs = 0.5, n = 2, zs = 0.004, 
                         Cp = 1005) {
  
  # constants
  Gsc <- 0.082
  k <- 0.41
  Rd <- 287.058
  Rv <- 461.495
  d <- 0
  
  # calculate environmental terms
  D <- RHtoVPD(RH = RH, TdegC = Ta, Pa = Pa)
  Tw <- wetbulb.temp(Tair = Ta, pressure = Pa * 1000, VPD = D)
  esTw <- 0.6108 * exp((17.27 * Tw) / (Tw + 238.3))
  ea <- esTw - y * Pa * (Ta - Tw)
  s <- (4098 * (0.6108 * exp((17.27 * Ta) / (Ta + 237.3)))) / (Ta + 237.3)^2 # Eq. 18
  ps <- 6.1078 * 10^((7.5 * Ta) / (Ta + 237.3))
  pv <- ps * RH
  pd <- (Pa * 1000) - pv
  p <- (pd / (Rd * (Ta + 273.15))) + (pv / (Rv * (Ta + 273.15)))
  
  # calculate distance from the earth to the sun, solar decimation, and sunset solar angle
  dr <- 1 + 0.033 * cos(J * ((2 * pi) / 365))
  sd <- 0.409 * sin(J * ((2 * pi) / 365) - 1.39)
  ws <- (pi/2) - atan((-tan(l) * tan(sd)) / (1 - ((tan(l))^2) * ((tan(sd))^2)^0.5))
  # calculate extraterrestrial radiation
  Ra <- (24 * 60 / pi) * Gsc * dr * 
    (ws * sin(l) * sin(sd) + cos(l) * cos(sd) * sin(ws))
  N <- ws * (24 / pi)
  # calculate shortwave radiation
  Rs <- ((as + bs * (ns / N)) * Ra)
  Rso <- (0.75 + (elev * 2e-5)) * Ra
  # calculate net shortwave radiation
  Rns <- (1 - alb) * Rs
  # calculate net longwave radiation
  Rnl <- (((Tmax + Tmin) / 2) * (0.34 - 0.14 * sqrt(ea)) * (1.35 * (Rs/Rso) - 0.35))
  # calculate net radiation
  Rn_MJ <- (Rns + Rnl)
  # convert from MJ m-2 day-1 to W m-2
  Rn <- 11.82033 * Rn_MJ
  
  # calculate ground radiation, Eq. 23
  G <- (Rn * 0.2)

  # calculate resistance terms
  rs <- exp(8.206 - 4.2555 * (SM / SMs))
  ra <- ((log10((x - d) / zs))^2) / ((k^2) * ux) # Eq. 25

  # calculate evapotranspiration using the Penman-Monteith equation
  ET <- (s * (Rn - G) + (p * Cp * D) / ra) / (s + y * (1 + rs / ra))
  
  # output
  return(data.frame(Rn = Rn, G = G, rs = rs, ra = ra, ET = ET))
}
