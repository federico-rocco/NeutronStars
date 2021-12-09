# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 12:13:14 2021

@author: cosmo
"""

#class Implicit:
    
import numpy as np
from scipy import interpolate
from utils import *



    
class Polytropic:
    
    
    def __init__(self, k, gamma):
        self.kind = "PressureEnergyDensityPolytropic"
        self.k = k*(cgs_geom_dictionary["cgs"]["lenght"]["m"]**2)*(cgs_geom_dictionary["cgs"]["energy"]["geom"]**(-2/3)) 
        self.gamma = gamma
        
    def DensityFromPressure(self, pressure):
        density = (pressure/self.k)**(1/self.gamma)
        return density
        
    def PressureFromDensity(self, density):
        pressure = self.k*density**self.gamma
        return pressure

class Piecewise:
 

#going from the core to the outside

    
    def __init__(self, core_gammas, central_density): #gammas = gamma1,gamma2,gamma3
        self.kind = "PressureEnergyDensityPiecewise"
        self.gammas = [core_gammas[0], core_gammas[1], core_gammas[2], 1.356, 0.622, 1.287, 1.584] #gammaSly3,sly2,sly1,sly0
        self.kappas = [0, 0, 0, 3.998e-8, 5.326e1, 1.061e-6, 6.801e-9] #k3,k2,k1,ksly3,ksly2,ksly1,ksly0
        self.rhos = [central_density, 1e15, 10**(14.7), 2.7e14, 2.627e12, 3.783e11, 2.440e7, 0] #rhocentral,3,2,trans,sy3,sly2,sly1
        self.rhos = [ x*cgs_geom_dictionary["cgs"]["density"]["geom"] for x in self.rhos]
        self.trans_pressure = 10**(34.384)*cgs_geom_dictionary["cgs"]["pressure"]["geom"]
        
    def BuildK(self):
        
        self.kappas = [ (self.kappas[i]*cgs_geom_dictionary["cgs"]["lenght"]["m"]**(3*self.gammas[i] -1)
                                           *cgs_geom_dictionary["cgs"]["time"]["geom"]**(-2)
                                           *cgs_geom_dictionary["cgs"]["mass"]["geom"]**(1-self.gammas[i]) )
                            for i in range(len(self.kappas))]     
        
        k1 = self.trans_pressure/(self.rhos[3]**self.gammas[2])
        k2 = k1*self.rhos[2]**(self.gammas[2]-self.gammas[1])
        k3 = k2*self.rhos[1]**(self.gammas[1]-self.gammas[0])
        self.kappas[2] = k1
        self.kappas[1] = k2
        self.kappas[0] = k3
    
    def BuildPressures(self):
        self.pressures = [ self.PressureFromDensity(density) for density in self.rhos]
    
    def PressureFromDensity(self, density):
        i=0

        while density <= self.rhos[i]:
            if density == 0:
                break
            i+=1

        return self.kappas[i-1]*density**self.gammas[i-1]
    
    def DensityFromPressure(self, pressure):
        i=0

        while pressure <= self.pressures[i]:
            if i==7:
                break
            i+=1
            
        return (pressure/self.kappas[i-1])**(1/self.gammas[i-1])

"""
6 -> gammaSLy0, kSly0, rhoSly1
5 -> gammaSLy1, kSly1, rhoSly2
4 -> gammaSLy2, kSly2, rhoSly3
3 -> gammaSLy3, kSly3, trans_density
2 -> gamma1, k1, rho2
1 -> gamma2, k2, rho3
0 -> gamma3, k3, rhoCentral
"""        



class Implicit:
    
    
    
    def __init__(self):
        self.kind = "PressureEnergyDensityImplicit"        
        self.e0 = ((m_n**4)*(C_SI**5)/((np.pi**2)*(HBAR_SI**3)))*cgs_geom_dictionary['si']['energy_density']['geom']
        self.DensityFromPressure, self.PressureFromDensity = self.InterpolateSolution()
        print("built eos")
        
        
        
    def ImplicitDensity(self, x):       
        return (self.e0/8)*((2*x**3 + x)*np.sqrt(1.0 + x**2) - np.arcsinh(x))

    def ImplicitPressure(self, x):
        return (self.e0/24)*((2*x**3 - 3*x)*np.sqrt(1.0 + x**2) + 3*np.arcsinh(x))
    
    
    def InterpolateSolution(self):
        x_range = np.linspace(0,1e3,int(100000))
        p_column, e_column = self.ImplicitPressure(x_range), self.ImplicitDensity(x_range)

        return CubicSpline(p_column, e_column), CubicSpline(e_column, p_column)

        

    
    
    
    
    
    
    
    
    
        
        
        
        