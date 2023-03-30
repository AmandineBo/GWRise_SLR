# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 14:48:43 2022

@author: abo132
"""

import numpy as np
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import math

"""
Based on the example of the Pioneer Valley, north-eastern Australia, QLD
Given in Morgan, L. K. and A. D. Werner (2016). "Comment on “Closed-form analytical solutions for assessing the consequences of sea-level rise on groundwater resources in sloping coastal aquifers”: paper published in Hydrogeology Journal (2015) 23:1399–1413, by R. Chesnaux." Hydrogeology journal 24(5): 1325-1328.

Werner & Gallagher (2006) Characterisation of sea-water intrusion in the Pioneer Valley, Australia using hydrochemistry and three-dim
"""
# aquifer parameters
k = 100 # hydraulic conductivity, m/d
z0 = 25 # (aquifer thickness) depth to the aquifer base measured from mean sea level, m
rhof = 1000 # density of freshwater, kg/m^3
rhos = 1025 # density of saltwater, kg/m^3
delt=(rhos-rhof)/rhof

W = 0.11/365.25 # net recharge from RNF, [L/T], accounting for infiltration, evapotranspiration and distributed pumping; mm to m; yr to d, m/d 
#W = 0.00005

# For an observation well with a water level at 2 m above MSL and 2000 m distant from coast
# # Numerical Answer FF
# hb=	1.24941
# xb=	1000

# Numerical Answer FH
# hb=	1.24568
# xb=	1000
# hb=	1.6159
# xb=	1500
hb=	2
xb=	2000


#%% Zone 1 - Inland of the interface
# lateral flow at observed location (xb, hb), m^2/d
qb = ((((hb+z0)**2 -((1+delt)*z0**2))*k) - (W*(xb)**2)) /(2*xb)
# flow toward the coast, m^2/d

q0 = qb+(W*xb)
# discharge to the sea, m^2/d

xn=q0/W
#Distance from the coast to groundwater divide (no flow inland)

#Mixed Convection Ratio - Vulnerability Measure, if M >= 1 density-driven processes dominate and there is a high SWI vulnerability:
M = k*delt*(1+delt)*z0**2/(W*xn**2)

xt = xn*(1-np.sqrt(1-M))
preSLR_Ratio =xt/xn

#%%
#At distance x and with DeltaZ (SLR)

x= 500
DeltaZ=1

#%% Under flux-controlled conditions:  z0 with z0 + ∆z
h = np.sqrt((((2*q0*x)-(W*x**2))/k)+((1+delt)*z0**2))-z0

hPrim = np.sqrt((((2*q0*x)-(W*x**2))/k)+((1+delt)*(z0+DeltaZ)**2))-(z0+DeltaZ)

Deltah = (hPrim -h)+DeltaZ

print(Deltah)


#Distance from the coast to groundwater divide (no flow inland) doesn't change Pre SLR and Post SLR under flux-controlled conditions because q0 doesn't change.

#Mixed Convection Ratio - Vulnerability Measure, if M >= 1 density-driven processes dominate and there is a high SWI vulnerability:
M_FF = k*delt*(1+delt)*(z0+DeltaZ)**2/(W*xn**2)

xt_FF_SLR = xn*(1-np.sqrt(1-M_FF))

#%% Under head-controlled conditions
#Under head-controlled conditions, z0 is replaced with z0 + ∆z and q0 is replaced with a new value determined from Equation 4 or 5.
qb_N = ((((hb+z0)**2 -((1+delt)*(z0+DeltaZ)**2))*k) - (W*(xb)**2)) /(2*xb)

q0_N = qb_N + (W*xb)

h = np.sqrt((((2*q0*x)-(W*x**2))/k)+((1+delt)*z0**2))-z0

xn_N=q0_N/W

print(xn-xn_N)

hPrim_N = np.sqrt((((2*q0_N*x)-(W*x**2))/k)+((1+delt)*(z0+DeltaZ)**2))-(z0+DeltaZ)


Deltah_N = (hPrim_N -h)+DeltaZ

print(Deltah_N)

#Mixed Convection Ratio - Vulnerability Measure, if M >= 1 density-driven processes dominate and there is a high SWI vulnerability:
M_FH = k*delt*(1+delt)*(z0+DeltaZ)**2/(W*xn_N**2)

xt_FH_SLR = xn_N*(1-np.sqrt(1-M_FH))

#%%
#Delta at the groundwater divide

hn = np.sqrt(((q0**2)/(W*k))+((rhos/rhof)*z0**2))-z0

hnPrim = np.sqrt(((q0**2)/(W*k))+((rhos/rhof)*(z0+DeltaZ)**2))-(z0+DeltaZ)

Deltahn = (hnPrim -hn)+DeltaZ

print(Deltahn)








