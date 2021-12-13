# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 17:10:23 2021

@author: cosmo
"""

#selezione eq stato
from project.eqs_state import Piecewise
from project.NSclass import NeutronStar as ns
from utils import cgs_geom_dictionary
import matplotlib.pyplot as plt
import numpy as np
import os


#creating results directory
script_dir = os.path.dirname(__file__)
output_dir = os.path.join(script_dir, "NSOutput/Read/")
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)
    


#creazione equazione di stato
cv = 3.5e15
chosen_state = Piecewise("SLy")
chosen_state.BuildK()
chosen_state.BuildPiecewise()
#chosen_state.BuildPressures()
print("Built equation of state")


#creating star
star = ns("ReadNS", chosen_state)




#Mass and pressure vs radius


#solving system

r_TOV,m_TOV,p_TOV = star.star_solver(star.TOV_eqs, cv, "density", "cgs")


iterations = r_TOV.size
print("TOV:",iterations, " iterations executed; mass = ", m_TOV[-1],"solar masses; total radius =", r_TOV[-1], "km")
R_TOV = r_TOV[-1]
M_TOV= m_TOV[-1]


#plot

fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.text(0.5, 1.07, "P(r) & m(r) of a Read neutron star", horizontalalignment='center', fontsize=12, transform = ax.transAxes)
ax.plot(r_TOV, p_TOV, color="black", linestyle="-", linewidth=2,  label = 'P TOV')
ax.set_xlabel('r [km]',fontsize=14)
ax.set_ylabel(r'P [$dyne/cm^2$]', fontsize=14)
ax.minorticks_on()

ax2 = ax.twinx()
ax2.minorticks_on()
ax2.plot(r_TOV, m_TOV, color="black", linestyle="-.", label = 'm TOV')
ax2.plot(R_TOV, M_TOV, marker = 'o', color='red', linestyle="", label='NS TOV mass')
ax2.set_ylabel(r"m [$M_{\odot}$]",fontsize=14)

fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(output_dir+'read_m&p-vs-radius.pdf', format='pdf', dpi=1000)


#sequence of 200 stars for Mass vs Radius and Mass/Radius vs central pressure 
 
d0_min = 13.5
d0_max = 16.5
densities = np.logspace(d0_min, d0_max, 200)
print("Solving 200 stars")
R_star_TOV, M_star_TOV = star.mass_vs_radius(densities, star.TOV_eqs, "density", "cgs")


# plot Mass vs Radius

fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.title("Mass-Radius of a Read NS")
ax.plot(R_star_TOV, M_star_TOV, color="black", linestyle="-.", linewidth=2,  label = 'TOV')
ax.set_xlabel('R [km]',fontsize=14)
ax.set_ylabel(r"M [$M_{\odot}$]", fontsize=14)
ax.minorticks_on()

fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(output_dir+'read_mass-vs-radius.pdf', format='pdf', dpi=1000)

# plot Mass/Radius vs central pressure 

fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.title("Mass/Radius vs Central Pressure in a Read NS")
ax.set_xscale('log')
ax.minorticks_on()
ax.plot(densities, M_star_TOV, color="black", linestyle="-.", linewidth=2,  label = 'M-TOV')
ax.set_xlabel('rho0 [$dyne/cm^2$]',fontsize=14)
ax.set_ylabel(r"M [$M_{\odot}$]", fontsize=14)

ax2 = ax.twinx()
ax2.set_xscale('log')
ax2.minorticks_on()
ax2.plot(densities, R_star_TOV, color="black", linestyle=":", linewidth=2,  label = 'R-TOV')
ax2.set_ylabel(r"R [$km$]",fontsize=14)

fig.legend(loc='upper center', bbox_to_anchor=(0.5,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(output_dir+'read_mr-vs-p0.pdf', format='pdf', dpi=1000)
