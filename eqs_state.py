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
    
    a = 0.0
    def __init__(self, k, gamma):
        self.kind = "PressureEnergyDensityPolytropic"
        self.k = k 
        self.gamma = gamma
        self.n = 1/(self.gamma-1)
        
    def DensityFromPressure(self, pressure):
        density = (pressure/self.k)**(1/self.gamma)
        return density
        
    def PressureFromDensity(self, density):
        #print("computing with", self.k, self.gamma)
        pressure = self.k*density**self.gamma
        return pressure
    
    def EdenFromPressure(self, pressure):
        density = self.DensityFromPressure(pressure)
        eden = self.EdenFromDensity(density)
        return eden

    def EdenFromDensity(self,density):
        eden = (1+self.a)*density + self.n*self.k*(density**self.gamma)
        return eden

class Piecewise:
 

#going from the core to the outside

    
    def __init__(self, key): 
        parameters = eos_library[key]
        self.kind = "PressureEnergyDensityPiecewise"
        self.gammas = [1.58425, 1.28733, 0.62223, 1.35692, parameters[1], parameters[2], parameters[3]] #gammaSly0,sly1,sly2,sly3,gamma1,2,3
        self.kappas = [6.80110e-9, 1.06186e-6, 5.32697e1, 3.99874e-8] #ksly0,sly1,sly2,sly3
        self.densities = [1e3,  2.44034e7, 3.78358e11, 2.62780e12, 2.7e14, 10**(14.7), 1e15] #rhosly1,2,3,trans,2,3,central        
        self.densities = [ x*cgs_geom_dictionary["cgs"]["density"]["geom"] for x in self.densities]
        self.trans_pressure = (10**parameters[0])*cgs_geom_dictionary["cgs"]["pressure"]["geom"] #k1*ptrans**Gamma1
        self.layers = []
        self.edens = []
        self.pressures = []
        
    def BuildK(self):
        
        self.kappas = [ (C_CGS**2)*(self.kappas[i]*cgs_geom_dictionary["cgs"]["lenght"]["m"]**(3*self.gammas[i] -1)
                                           *cgs_geom_dictionary["cgs"]["time"]["geom"]**(-2)
                                           *cgs_geom_dictionary["cgs"]["mass"]["geom"]**(1-self.gammas[i]) )
                            for i in range(len(self.kappas))]     
        
        k1 = self.trans_pressure/(self.densities[4]**self.gammas[4])
        k2 = k1*self.densities[5]**(self.gammas[4]-self.gammas[5])
        k3 = k2*self.densities[6]**(self.gammas[5]-self.gammas[6])
        self.kappas.append(k1)
        self.kappas.append(k2)
        self.kappas.append(k3)
        print("kappas:",self.kappas,"gammas:",self.gammas)
        #self.kappas = [1.0770248668221202e+17, 1.9183640140950262e+18, 2.7747766144672896e+18, 16.373064274559898, 4.0848760905320256e-08, 9.101658781649894, 851304.850457802]
    
    def BuildPressures(self):
        self.pressures = [ self.PressureFromDensity(density) for density in self.densities]
        #self.pressures[3] = 1.0686386273362236e-13
        print("densities:",self.densities,"pressures:",self.pressures)
        
    def BuildPiecewise(self):
        for k, gamma in zip(self.kappas, self.gammas):
            
            self.layers.append(Polytropic(k, gamma))
        self.densities[4] = (self.kappas[4]/self.kappas[3])**(1/(self.gammas[3]-self.gammas[4]))
            #print(k,gamma,self.layers[-1])
        
        prev_trope=None
        for (layer, density) in zip(self.layers, self.densities):

                if not( prev_trope == None ):
                    layer.a = self.set_trans_const( prev_trope, layer, density )
                else:
                    density = 0.0

                eden = layer.EdenFromDensity(density)
                pressure = layer.PressureFromDensity(density)

                # eden and pressures corresponding to the densities transition values

                self.pressures.append( pressure )
                self.edens.append( eden )

                prev_ed = eden
                prev_tr = density
                prev_trope = layer

    def set_trans_const(self, prev_trope, trope, transition):
        
        return prev_trope.a + (prev_trope.k/(prev_trope.gamma - 1))*transition**(prev_trope.gamma-1) - (trope.k/(trope.gamma - 1))*transition**(trope.gamma-1)
        
        
        
        
        
    def PressureFromDensity(self, density):        
        layer = self.FindLayer(density, "density")
        #print("for density",density,"found",layer.PressureFromDensity(density) )
        return layer.PressureFromDensity(density) 


    def DensityFromPressure(self, pressure,index="False"):                
        layer = self.FindLayer(pressure, "pressure")
        return layer.DensityFromPressure(pressure)     

    def EdenFromPressure(self, pressure,index="False"):                
        layer = self.FindLayer(pressure, "pressure")
        return layer.EdenFromPressure(pressure)    
    def FindLayer(self, value, value_type):
        #print("search press")
        if value_type == "pressure":
            #print("searching")
            if value >= self.pressures[-1]:
                return self.layers[-1]
            for x in range( len(self.pressures) - 1):
                if self.pressures[x] <= value < self.pressures[x+1]:
                    #print("found pressure", x)
                    return self.layers[x]
                
        elif value_type == "density":
            #print("searchingd")
            if value >= self.densities[-1]:
                return self.layers[-1]
            for x in range( len(self.densities) - 1):
                if self.densities[x] <= value < self.densities[x+1]:
                    #print("found dens", x)
                    return self.layers[x] 
                
                
                
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

        

    
    
    
    
    
    
    
    
    
        
        
        
        