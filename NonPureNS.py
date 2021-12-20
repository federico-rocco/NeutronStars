# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 17:10:23 2021

@author: cosmo
"""

#selezione eq stato
from eqs_state import Piecewise
from NSclass import NeutronStar as ns
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
chosen_state.build_piecewise()
print("Built equation of state")


#creating star
star = ns("ReadNS", chosen_state)




#Mass and pressure vs radius


#solving system

r_tov,m_tov,p_tov = star.star_solver(star.tov_eqs, cv, "density", "cgs")


iterations = r_tov.size
print("tov:",iterations, " iterations executed; mass = ", m_tov[-1],"solar masses; total radius =", r_tov[-1], "km")
R_tov = r_tov[-1]
M_tov= m_tov[-1]


#plot

fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.text(0.5, 1.07, "P(r) & m(r) of a Read neutron star", horizontalalignment='center', fontsize=12, transform = ax.transAxes)
ax.plot(r_tov, p_tov, color="black", linestyle="-", linewidth=2,  label = 'P tov')
ax.set_xlabel('r [km]',fontsize=14)
ax.set_ylabel(r'P [$dyne/cm^2$]', fontsize=14)
ax.minorticks_on()

ax2 = ax.twinx()
ax2.minorticks_on()
ax2.plot(r_tov, m_tov, color="black", linestyle="-.", label = 'm tov')
ax2.plot(R_tov, M_tov, marker = 'o', color='red', linestyle="", label='NS tov mass')
ax2.set_ylabel(r"m [$M_{\odot}$]",fontsize=14)

fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(output_dir+'read_m&p-vs-radius.pdf', format='pdf', dpi=1000)


#sequence of 200 stars for Mass vs Radius and Mass/Radius vs central pressure 
 
d0_min = 14
d0_max = 16.5
densities = np.logspace(d0_min, d0_max, 20)
print("Solving 200 stars")
R_star_tov, M_star_tov = star.mass_vs_radius(densities, star.tov_eqs, "density", "cgs")

# plot Mass vs Radius

fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.title("Mass-Radius of a Read NS")
ax.plot(R_star_tov, M_star_tov, color="black", linestyle="-.", linewidth=2,  label = 'tov')
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
ax.plot(densities, M_star_tov, color="black", linestyle="-.", linewidth=2,  label = 'M-tov')
ax.set_xlabel('$p_0$ [$dyne/cm^2$]',fontsize=14)
ax.set_ylabel(r"M [$M_{\odot}$]", fontsize=14)

ax2 = ax.twinx()
ax2.set_xscale('log')
ax2.minorticks_on()
ax2.plot(densities, R_star_tov, color="black", linestyle=":", linewidth=2,  label = 'R-tov')
ax2.set_ylabel(r"R [$km$]",fontsize=14)

fig.legend(loc='upper center', bbox_to_anchor=(0.5,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(output_dir+'read_mr-vs-p0.pdf', format='pdf', dpi=1000)