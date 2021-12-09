# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 19:02:16 2021

@author: cosmo
"""

import numpy as np
from project.utils import ODEsolver, cgs_geom_dictionary


class NeutronStar:
    
    def __init__(self, NS_type, eq_state):
        
        self.kind = NS_type
        self.eos = eq_state
        
    
    def TOV_eqs(self, r, y):

        m,p=y
        eden = self.eos.DensityFromPressure(p)
        dm = 4*np.pi*eden*r**2
        dp = - (eden+p)*(m + 4*np.pi*r**3*p)/(r*(r-2*m))
        dy = [dm, dp]
        
        return dy

    def Newton_eqs(self, r, y):

        m,p=y
        eden = self.eos.DensityFromPressure(p)
        dm = 4*np.pi*eden*r**2
        dp = - (eden*m)/(r**2)
        dy = [dm, dp]

        return dy
    
    def star_solver(self, eq_type, central_value, value_type="pressure", unit_type="geom"):

        if value_type == "density":
            central_value = self.eos.PressureFromDensity(central_value)
        
        if unit_type == "cgs":
            central_value = central_value*cgs_geom_dictionary["cgs"]["pressure"]["geom"]
        elif unit_type == "si":
            central_value = central_value*cgs_geom_dictionary["si"]["pressure"]["geom"]    


        solutions = ODEsolver(eq_type, central_value, self.eos)
        
        r_out = np.array([])
        m_out = np.array([])
        p_out = np.array([])

        for i in range(solutions[0].size):
            r_out = np.append(r_out, solutions[0][i]*cgs_geom_dictionary["geom"]["lenght"]["km"])
            m_out = np.append(m_out, solutions[1][i]*cgs_geom_dictionary["geom"]["mass"]["m_sol"])
            p_out = np.append(p_out, solutions[2][i]*cgs_geom_dictionary["geom"]["pressure"]["cgs"])

        return r_out, m_out, p_out
    
    
    def mass_vs_radius(self, pressures, eq_type):
        
        radii = np.array([])
        masses = np.array([])
        
        import tqdm
        import progressbar
        from time import sleep
        bar = progressbar.ProgressBar(maxval=20, \
            widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        #bar.start()

        
        for i in tqdm.tqdm(range(pressures.size)):
        #for i in range(pressures.size):
            #bar.update(i)
            solutions = self.star_solver(eq_type, pressures[i], value_type="pressure", unit_type="cgs")
            radii = np.append(radii, solutions[0][-1])
            masses = np.append(masses, solutions[1][-1])
        #bar.finish()
        
        return radii, masses
        
            
            