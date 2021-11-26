# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 15:47:45 2021

@author: cosmo
"""

#selezione eq stato
from project.eqs_state import eqs_state
eqs_state.select()

#soluzione sistema
from solve import solve
radius, mass, pressure = solve()

#plot
import matplotlib.pyplot as plt
plt.plot(radius,mass,radius,pressure)
plt.show()