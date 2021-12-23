# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 12:19:57 2021

@author: cosmo
"""

#selezione eq stato
from eqs_state import Implicit
from NSclass import NeutronStar as ns
import matplotlib.pyplot as plt
import numpy as np
import os


#creating results directory
script_dir = os.path.dirname(__file__)
output_dir = os.path.join(script_dir, "NSOutput/Rel/")
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)
    
    
#creating equation of state
chosen_state = Implicit()
print("eos made")

#creating star
star = ns("RelPureNS", chosen_state)

#Mass and pressure vs radius


#solving system

cv = 3.5*10**35
r_newton, m_newton, p_newton = star.star_solver(star.newton_eqs, cv, "pressure", "cgs")
r_tov, m_tov, p_tov = star.star_solver(star.tov_eqs, cv, "pressure", "cgs")

iterations = r_newton.size
print("Newton:", iterations, " iterations executed; mass = ", m_newton[-1],"solar masses; total radius =", r_newton[-1], "km")
R_newton = r_newton[-1]
M_newton= m_newton[-1]

iterations = r_tov.size
print("tov:",iterations, " iterations executed; mass = ", m_tov[-1],"solar masses; total radius =", r_tov[-1], "km")
R_tov = r_tov[-1]
M_tov= m_tov[-1]


#plot

fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.text(0.5, 1.07, "P(r) & m(r) of a pure relativistic neutron star", horizontalalignment='center', fontsize=12, transform = ax.transAxes)
ax.plot(r_newton, p_newton, color="blue", linestyle="-", linewidth=1, label = 'P Newton')
ax.plot(r_tov, p_tov, color="black", linestyle="-", linewidth=2,  label = 'P tov')
ax.set_xlabel('r [km]',fontsize=14)
ax.set_ylabel(r'P [$dyne/cm^2$]', fontsize=14)
ax.minorticks_on()

ax2 = ax.twinx()
ax2.minorticks_on()
ax2.plot(r_newton, m_newton,color="blue", linestyle=":", label = 'm Newton')
ax2.plot(r_tov, m_tov, color="black", linestyle="-.", label = 'm tov')
ax2.plot(R_newton, M_newton, marker = 'o', linestyle="", color='green', label='NS Newton mass')
ax2.plot(R_tov, M_tov, marker = 'o', color='red', linestyle="", label='NS tov mass')
ax2.set_ylabel(r"m [$M_{\odot}$]",fontsize=14)

fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(output_dir+'relns_m&p-vs-radius.pdf', format='pdf', dpi=1000)


#sequence of 200 stars for Mass vs Radius and Mass/Radius vs central pressure 

p0_min = 3.5e34
p0_max = 3.5e38
pressures = np.linspace(p0_min, p0_max, 1000)
R_star_tov = np.zeros(pressures.size)
M_star_tov = np.zeros(pressures.size)
R_star_newton = np.zeros(pressures.size)
M_star_newton = np.zeros(pressures.size)
print("Solving 1000 stars")
if __name__ == '__main__':
    R_star_tov, M_star_tov = star.mass_vs_radius(pressures, star.tov_eqs, "pressure", "cgs")
    R_star_newton, M_star_newton = star.mass_vs_radius(pressures, star.newton_eqs, "pressure", "cgs")


# plot Mass vs Radius

fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.title("Mass-Radius of a pure relativistic NS")
ax.plot(R_star_newton, M_star_newton, color="blue", linestyle="-.", linewidth=1, label = 'Newton')
ax.plot(R_star_tov, M_star_tov, color="black", linestyle="-.", linewidth=2,  label = 'tov')
ax.set_xlabel('R [km]',fontsize=14)
ax.set_ylabel(r"M [$M_{\odot}$]", fontsize=14)
ax.minorticks_on()

fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(output_dir+'relns_mass-vs-radius.pdf', format='pdf', dpi=1000)

# plot Mass/Radius vs central pressure 

fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.title("Mass/Radius vs Central Pressure in a pure relativistic NS")
ax.set_xscale('log')
ax.minorticks_on()
ax.plot(pressures, M_star_newton, color="blue", linestyle="-.", linewidth=1, label = 'M-Newton')
ax.plot(pressures, M_star_tov, color="black", linestyle="-.", linewidth=2,  label = 'M-tov')
ax.set_xlabel('$p_0$ [$dyne/cm^2$]',fontsize=14)
ax.set_ylabel(r"M [$M_{\odot}$]", fontsize=14)

ax2 = ax.twinx()
ax2.set_xscale('log')
ax2.minorticks_on()
ax2.plot(pressures, R_star_newton, color="blue", linestyle=":", linewidth=1, label = 'R-Newton')
ax2.plot(pressures, R_star_tov, color="black", linestyle=":", linewidth=2,  label = 'R-tov')
ax2.set_ylabel(r"R [$km$]",fontsize=14)

fig.legend(loc='upper center', bbox_to_anchor=(0.5,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(output_dir+'relns_mr-vs-p0.pdf', format='pdf', dpi=1000)