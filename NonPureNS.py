# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 17:10:23 2021

@author: cosmo
"""

#selezione eq stato
from project.eqs_state import Piecewise
from project.NSclass import NeutronStar as ns
from utils import cgs_geom_dictionary
import numpy as np



#creazione equazione di stato
chosen_state = Piecewise([2.159,  2.189,  2.227], trans_pressure=10**34.38, presCentral=1e40)
chosen_state.BuildK()
print("eos made")

#creazione stella
star = ns("NonPureNS", chosen_state)

#soluzione sistema
cv = 3.5e15*cgs_geom_dictionary["cgs"]["density"]["geom"]
r_TOV,m_TOV,p_TOV = star.star_solver(star.TOV_eqs, cv, "density", "geom")

iterations = r_TOV.size
print("TOV:",iterations, " iterations executed; mass = ", m_TOV[-1],"solar masses; total radius =", r_TOV[-1], "km")
R_TOV = r_TOV[-1]
M_TOV= m_TOV[-1]


#plot
import matplotlib.pyplot as plt


fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.text(0.5, 1.07, "P(r) & m(r) of a non pure relativistic neutron star",
         horizontalalignment='center',
         fontsize=12,
         transform = ax.transAxes)

ax.plot(r_TOV, p_TOV, color="black", linestyle="-", linewidth=2,  label = 'P TOV')
ax.set_xlabel('r [km]',fontsize=14)
ax.set_ylabel(r'P [$dyne/cm^2$]', fontsize=14)
ax.minorticks_on()

ax2 = ax.twinx()
ax2.minorticks_on()

ax2.plot(r_TOV, m_TOV, color="black", linestyle="-.", label = 'm TOV')

ax2.plot(R_TOV, M_TOV, marker = 'o', color='red', linestyle="", label='NS TOV mass')
ax2.set_ylabel(r"m [$M_{\odot}$]",fontsize=14)
