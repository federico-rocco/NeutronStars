# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 12:13:14 2021

@author: cosmo
"""

#class Implicit:
    
import numpy as np    
    
class Polytropic:
    
    
    def __init__(self, k, gamma):
        self.kind = "PressureEnergyDensityPolytropic"
        self.k = k
        self.gamma = gamma
        
    def DensityFromPressure(self, pressure):
        density = (pressure/self.k)**(1/self.gamma)
        return density
        
    def PressureFromDensity(self, density):
        pressure = self.k*density**self.gamma
        return pressure

class Piecewise:
 

#going from the core to the outside

    
    def __init__(self, trans_density, core_gammas, rhoCentral = "bho"): #gammas = gamma1,gamma2,gamma3
        self.kind = "PressureEnergyDensityPiecewise"
        self.gammas = [core_gammas, 1.356, 0.622, 1.287, 1.584] #gammaSly3,sly2,sly1,sly0
        self.rhos = [rhoCentral, 10**(5), 10**(4.7), trans_density, 2.627*10**12, 3.783*10**11, 2.440*10**7] #rhocentral,3,2,trans,sy3,sly2,sly1
        self.kappas = [3.998*10**(-8), 5.326*10**(1), 1.061*10**(-6), 6.801*10**(-9)] #ksly3,sly2,sly1,sly0
        
    def BuildK(self):
        for i in range(2, 0, -1):
            new_k = self.PressureFromDensity(self.rhos[i+1])/(self.rhos[i+1]**self.gammas[i])
            self.kappas.insert(0,new_k)
    
    def PressureFromDensity(self, density):
        i=0
        while density < self.rhos[i+1]:
            i+=1
        return self.kappas[i]*density**self.gammas[i]

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
        self.InterpolatedPressure, self.InterpolatedDensity = self.InterpolateSolution()
        
        
    def density(self, x):
        return (1/8)*((2*x**3 + x)*np.sqrt(1 + x**2) - np.arcsinh(x))

    def pressure(self, x):
        return (1/24)*((2*x**3 - 3*x)*np.sqrt(1 + x**2) + 3*np.arcsinh(x))
    
    
    def InterpolateSolution(self):
        a=2
        b=3
        from solve import bisection
        eden_guess = np.linspace(0, b-a, b-a)
        
        k_fermi = []
        for i in range(b-a):
            k_fermi.append(bisection(self.density, a, b, eden_guess[i]))

        p = []
        for j in range(b-a):
            p.append(self.pressure(k_fermi[j]))

        return np.interpolate.interp1d(eden_guess, p, kind = 3), np.interpolate.interp1d(p, eden_guess, kind = 3)
    
    def DensityFromPressure(self, pressure):
        return self.InterpolatedDensity(pressure)
  
    def PressureFromDensity(self, density):
        return self.InterpolatedPressure(density)    
    
    
    
    
    
    
    
    
    
        
        
        
        