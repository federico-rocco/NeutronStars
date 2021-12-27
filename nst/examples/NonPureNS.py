# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 17:10:23 2021

@author: cosmo
"""

import matplotlib.pyplot as plt
import numpy as np
import os

from nst.eqs_state import Piecewise
from nst.NSclass import NeutronStar as ns



#creating results directory
script_dir = os.path.dirname(__file__)
output_dir = os.path.join(script_dir, "NSOutput/Read/")
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)
    


#creazione equazione di stato
central = 3.5e15
chosen_state = Piecewise("SLy")
chosen_state.build_piecewise()

#creating star
star = ns("ReadNS", chosen_state)



#Mass and pressure vs radius

#solving system
r_tov,m_tov,p_tov = star.star_solver(star.tov_eqs, central, "density", "cgs")

iterations_tov = r_tov.size
R_tov = r_tov[-1]
M_tov= m_tov[-1]

print("=========================================================")
print("Solved SLy Neutron Star.")
print("Central density: ", central, " g/cm^3")
print("---------------------------------------------------------")
print("TOV:",iterations_tov, " iterations executed; mass = ", m_tov[-1],"solar masses; total radius =", r_tov[-1], "km")
print("=========================================================")



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



#sequence of 1000 stars for Mass vs Radius and Mass/Radius vs central pressure  
d_min = 13.9
d_max = 16.5
densities = np.logspace(d_min, d_max, 500)
print("Simulating a sequence of 500 Read Neutron Stars with densities between ", d_min, " and ", d_max, "g/cm^3")
print("---------------------------------------------------------")

print("Using TOV equations")
R_star_tov = np.zeros(densities.size)
M_star_tov = np.zeros(densities.size)
if __name__ == '__main__':
    R_star_tov, M_star_tov = star.mass_vs_radius(densities, star.tov_eqs, value_type="density", unit_type="cgs")



#plot Mass vs Radius
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



#plot Mass/Radius vs central density 
fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.title("Mass/Radius vs central density in a Read NS")
ax.set_xscale('log')
ax.minorticks_on()
ax.plot(densities, M_star_tov, color="black", linestyle="-.", linewidth=2,  label = 'M-tov')
ax.set_xlabel('$d_0$ [$g/cm^3$]',fontsize=14)
ax.set_ylabel(r"M [$M_{\odot}$]", fontsize=14)

ax2 = ax.twinx()
ax2.set_xscale('log')
ax2.minorticks_on()
ax2.plot(densities, R_star_tov, color="black", linestyle=":", linewidth=2,  label = 'R-tov')
ax2.set_ylabel(r"R [$km$]",fontsize=14)

fig.legend(loc='upper center', bbox_to_anchor=(0.5,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(output_dir+'read_mr-vs-d0.pdf', format='pdf', dpi=1000)



print("=========================================================")
print("Plots of the single star: Mass and Pressure vs Radius") 
print("Plots of the sequence: Mass vs Radius and Mass/Radius vs central density")
print("available in ", output_dir)
print("=========================================================")