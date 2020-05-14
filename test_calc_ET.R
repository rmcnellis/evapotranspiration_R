# test dual-source ET predicitons

source('calc_ET.R')

test_ET <- calc_ET(Pa = 90, Ta = 25, RH = 40, J = 180, l = 0.58, elev = 900, ns = 6, 
                   alb = 0.2, Tmax = 30, Tmin = 20, L = 2, h = 0.3, ux = 2, x = 2, gl = 400, 
                   SM = 0.1)

View(test_ET)

test_ETs <- calc_ETs(Pa = 90, Ta = 25, RH = 40, J = 180, l = 0.58, elev = 900, ns = 6, 
                     alb = 0.2, Tmax = 30, Tmin = 20, ux = 2, SM = 0.1)

View(test_ETs)
