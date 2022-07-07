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

& Werner AD, Gallagher MR (2006) Characterisation of sea-water intrusion in the Pioneer Valley, Australia using hydrochemistry and three-dim
"""
# aquifer parameters
k = 100 # hydraulic conductivity, m/d
z0 = 25 # (aquifer thickness) depth to the aquifer base measured from mean sea level, m
rhof = 1000 # density of freshwater, kg/m^3
rhos = 1025 # density of saltwater, kg/m^3
delt=(rhos-rhof)/rhof

W = 0.11/365.25 # net recharge from RNF, [L/T], accounting for infiltration, evapotranspiration and distributed pumping; mm to m; yr to d, m/d 

# For c well 2 m above MSL and 2000 m distant from coast
hb=	2
xb=	2000

#%%
# lateral flow at observed location (xb, hb), m^2/d
qb = ((((hb+z0)**2 -((1+delt)*z0**2))*k) - (W*(xb)**2)) /(2*xb)
# flow toward the coast, m^2/d

q0 = qb+(W*xb)
# discharge to the sea, m^2/d

xn=q0/W
#Distance from the coast to groundwater divide (no flow inland)

#%%
#At distance x and with DeltaZ (SLR)

x= 1000
DeltaZ=1

#%%
h = np.sqrt((((2*qb*x)+(W*x**2))/k)+((1+delt)*z0**2))-z0

hPrim = np.sqrt((((2*qb*x)+(W*x**2))/k)+((1+delt)*(z0+DeltaZ)**2))-(z0+DeltaZ)

Deltah = (hPrim -h)+DeltaZ

print(Deltah)


#%%
#Delta at the groundwater divide

hn = np.sqrt(((q0**2)/(W*k))+((rhos/rhof)*z0**2))-z0

hnPrim = np.sqrt(((q0**2)/(W*k))+((rhos/rhof)*(z0+DeltaZ)**2))-(z0+DeltaZ)

Deltahn = (hnPrim -hn)+DeltaZ

print(Deltahn)








