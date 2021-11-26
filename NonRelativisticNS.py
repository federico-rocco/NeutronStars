# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 15:47:45 2021

@author: cosmo
"""

#selezione eq stato
from project.eqs_state import Polytropic
from project.NSclass import NeutronStar as ns
#from project.solve import *



#creazione equazione di stato
k=1
gamma=1
chosen_state = Polytropic(k, gamma)

#creazione stella
star = ns("NonrelPureNS", 0, chosen_state)

#soluzione sistema
cv = 1000.0
kind='Newton'
r_newt, m_newt, p_nwet = star.star_solver(star.Newton_eqs, cv)
r_tov, m_tov, p_tov = star.star_solver(star.TOV_eqs, cv)

#plot
"""import matplotlib.pyplot as plt
plt.plot(radius,mass,radius,pressure)
plt.show()"""