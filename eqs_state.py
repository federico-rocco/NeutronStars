# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 12:13:14 2021

@author: cosmo
"""

#class Implicit:
    
    
    
class Polytropic:
    
    
    def __init__(self, k, gamma):
        self.kind = "PressureEnergyDensityPolytropic"
        self.k = k
        self.gamma = gamma
        
    def energy_density(self, pres):
        en_den = (pres/self.k)**(1/self.gamma)
        return en_den
        

class Piecewise:
 


#going from the core to the outside

    
    def __init__(self, trans_density, core_gammas, rhoCentral = "bho"): #gammas = gamma1,gamma2,gamma3
        self.kind = "PressureEnergyDensityPiecewise"
        self.gammas = [core_gammas, gammaSly3, gammaSLy2, gammaSly1, gammaSly0]
        self.rhos = [rhoCentral, rho3, rho2, trans_density, rhoSly3, rhoSly2, rhoSly1]
        self.kappas = [kSly3, kSly2, kSly1, kSly0]
        
    def BuildK(self):
        for i in range(2, 0, -1):
            new_k = self.PressureFromDensity(self.rhos[i+1])/(self.rhos[i+1]**self.gammas[i])
            self.kappas.insert(0,new_k)
        
    def PressureFromDensity(self, density):
        i=-1
        while density < self.rhos[i]:
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
        
        
        
        
        
        