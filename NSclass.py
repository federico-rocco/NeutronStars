# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 19:02:16 2021

@author: cosmo
"""

import numpy as np
from project.utils import solve, cgs_geom_dictionary
from scipy.integrate import solve_ivp

def CGS_converter(a):
    return a

class NeutronStar:
    
    def __init__(self, NS_type, eq_state):
        
        self.kind = NS_type
        self.eos = eq_state
        
    
    def TOV_eqs(self, r, y):

        m,p=y
        #print(p)
        eden = self.eos.DensityFromPressure(p)
        #print("eden=",eden,"m=",m,"p=",p)
        dm = 4*np.pi*eden*r**2
        dp = - (eden+p)*(m + 4*np.pi*r**3*p)/(r*(r-2*m))
        #print("dp=",dp)
        dy = [dm, dp]
        return dy

    def Newton_eqs(self, r, y):

        m,p=y
        eden = self.eos.DensityFromPressure(p)
        #print(eden)
        dm = 4*np.pi*eden*r**2
        dp = - (eden*m)/(r**2)
        dy = [dm, dp]
        #print("r=",r,"dm=",dm,"dp=",dp)
        if type(eden) == int or type(eden) == float:
            a=1
        else:
            print(type(eden),eden)
        return dy
    
    def star_solver(self, eq_type, central_value, value_type="pressure"):

        if value_type == "density":
            central_value = self.eos.PressureFromDensity(central_value)
            
        central_value = central_value*cgs_geom_dictionary["cgs"]["pressure"]["geom"]
        
        print("starting solving")
        solutions = solve(eq_type, central_value, self.eos)
        #solutions = solve_ivp(eq_type, [10**(-5), 10**6], [0, central_value])
        #print(len(central_value))
        #print(solutions)
        #solutions = solve_ivp(eq_type, (0.0001,10**40),(0,central_value))    
        r_out = np.array([])
        m_out = np.array([])
        p_out = np.array([])
        print(solutions[1])
        for i in range(solutions[0].size):
            r_out = np.append(r_out, solutions[0][i]*cgs_geom_dictionary["geom"]["lenght"]["km"])
            m_out = np.append(m_out, solutions[1][i]*cgs_geom_dictionary["geom"]["mass"]["m_sol"])
            p_out = np.append(p_out, solutions[2][i]*cgs_geom_dictionary["geom"]["pressure"]["cgs"])
            

        """r_out = solutions.t*cgs_geom_dictionary['geom']['lenght']['km'] 
        m_out = solutions.y[0]*cgs_geom_dictionary['geom']['mass']['m_sol']  
        p_out = solutions.y[1] *cgs_geom_dictionary['geom']['pressure']['cgs']"""
        
        return r_out, m_out, p_out
    
    
        

        
        
        
        
            
            