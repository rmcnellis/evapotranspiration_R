# test dual-source ET predicitons

source('calc_ET.R')

test_dual_ET <- calc_dual_ET(Ta = seq(5, 35, 5))

View(test_ET)

test_soil_ET <- calc_soil_ET(Ta = seq(5, 35, 5))

View(test_ETs)
