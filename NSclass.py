# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 19:02:16 2021

@author: cosmo
"""

import numpy as np
import project.eqs_state
from project.solve import solve

def CGS_converter(a):
    return a

class NeutronStar:
    
    def __init__(self, NS_type, value_type, eq_state):
        
        self.kind = NS_type
        self.eos = eq_state
        
    #qualcosa col valuetype
    
    def TOV_eqs(self, r, m, p):

        eden = self.eos.energy_density(p)
        dm = 4*np.pi*eden*r**2
        dp = - (eden+p)*(m + 4*np.pi*r**3*p)/(r*(r-2*m))
        dy = [dm, dp]
        return dy

    def Newton_eqs(self, r, m, p):

        eden = self.eos.energy_density(p)
        dm = 4*np.pi*eden*r**2
        dp = - (eden*m)/(r**2)
        dy = [dm, dp]
        return dy
    
    def star_solver(self, eq_type, central_value):

        solutions = solve(eq_type, central_value, self.eos)
            
        r_out = CGS_converter(solutions[0])
        m_out = CGS_converter(solutions[1])
        p_out = CGS_converter(solutions[2])
        
        return r_out, m_out, p_out
    
    
        
        
        
        
        
        
            
            