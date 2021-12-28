# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 12:13:14 2021

@author: cosmo
"""


import numpy as np
import nst.utils as ut
from nst.utils import cgs_geom_dictionary, eos_library, cubic_spline


   
class Polytropic:
    
    
    
    def __init__(self, k, gamma):
    
        """
        
        Class that builds a polytropic equation of state relating pressure and energy density of the kind p = k*epsilon^gamma

        Parameters
        ----------
        k : float
            tropic coefficient
        gamma : float
            adiabatic index
        n : float
            polytropic index
        a : float
            transition continuity constant 
        -------

        """
        
        self.kind = "PressureEnergyDensityPolytropic"
        self.k = k 
        self.gamma = gamma
        self.n = 1/(self.gamma-1)
        self.a = 0.0
        
        
        
    def density_from_pressure(self, pressure):
    
        """
        
        Method that returns the density given the pressure, as d = (p/k)^(1/gamma)

        Parameters
        ----------
        pressure : float
            input pressure

        Returns density : float
            computed density
        -------

        """
        
        density = (pressure/self.k)**(1/self.gamma)
        return density
        
    
    
    def pressure_from_density(self, density):
        
        """
        
        Method that returns the pressure given the density, as p = k*d^gamma

        Parameters
        ----------
        density : float
            input density

        Returns pressure : float
            computed pressure
        -------

        """
        
        pressure = self.k*density**self.gamma
        return pressure
    
    
    
    def eden_from_pressure(self, pressure):
        
        """
        
        Method that returns the energy density given the pressure, by calling density_from_pressure and eden_from_density

        Parameters
        ----------
        pressure : float
            input pressure

        Returns eden : float
            computed energy density
        -------

        """
        
        density = self.density_from_pressure(pressure)
        eden = self.eden_from_density(density)
        return eden



    def eden_from_density(self,density):
        
        """
        
        Method that returns the energy density given the density, as epsilon = (1+a)*d + n*k*d^gamma

        Parameters
        ----------
        density : float
            input density

        Returns eden : float
            computed energy density
        -------

        """

        eden = (1 + self.a)*density + self.n*self.k*(density**self.gamma)
        return eden



class Piecewise:
    
    
    
    def __init__(self, key):
    
        """
        
        Class that builds a piecewise equation of state relating pressure and energy density as p = k_i*epsilon^gamma_i, 
        where i represent the layer of the star. Each layer is a Polytropic object

        Parameters
        ----------
        key : string
            allows to choose a set of values (transition pressure, gamma_1, gamma_2, gamma_3) from a set of different models
            provided in eos_library in utils
        kappas : array of float
            tropic coefficients
        gammas : array of float
            adiabatic indexes
        trans_pressure : float
            pressure of transition between core and crust, needed to compute k1 in build_k
        densities : array of float
            densities of transition between layers
        pressures : array of float
            pressures of transition between layers, computed in build_piecewise
        edens : array of float
            energy densities of transition between layers, computed in build_piecewise
        layers : array of Polytropic objects
            polytropic equation of states for each layer
        
        -------
        Note: the piecewise eos is not complete with the constructor, build_piecewise needs to be called

        """
        
        parameters = eos_library[key]
        self.kind = "PressureEnergyDensityPiecewise"
        self.kappas = [6.80110e-9, 1.06186e-6, 5.32697e1, 3.99874e-8] #ksly0,sly1,sly2,sly3
        self.gammas = [1.58425, 1.28733, 0.62223, 1.35692, parameters[1], parameters[2], parameters[3]] #gammaSly0,sly1,sly2,sly3,gamma1,2,3
        self.trans_pressure = (10**parameters[0])*cgs_geom_dictionary["cgs"]["pressure"]["geom"] #k1*ptrans**Gamma1
        self.densities = [0,  2.44034e7, 3.78358e11, 2.62780e12, 2.7e14, 10**(14.7), 1e15] #rhosly1,2,3,trans,2,3,central        
        self.densities = [ x*cgs_geom_dictionary["cgs"]["density"]["geom"] for x in self.densities]
        self.pressures = []
        self.edens = []
        self.layers = []
   
        
   
    def build_k(self):
    
        """
        
        Method that computed the tropic coefficients of the core starting from the transition pressure
        -------

        """    
        
        self.kappas = [ (ut.C_CGS**2)*(self.kappas[i]*cgs_geom_dictionary["cgs"]["lenght"]["m"]**(3*self.gammas[i] -1)
                                           *cgs_geom_dictionary["cgs"]["time"]["geom"]**(-2)
                                           *cgs_geom_dictionary["cgs"]["mass"]["geom"]**(1-self.gammas[i]) )
                            for i in range(len(self.kappas))]     
        
        k1 = self.trans_pressure/(self.densities[4]**self.gammas[4])
        k2 = k1*self.densities[5]**(self.gammas[4]-self.gammas[5])
        k3 = k2*self.densities[6]**(self.gammas[5]-self.gammas[6])
        self.kappas.append(k1)
        self.kappas.append(k2)
        self.kappas.append(k3)

        

    def build_piecewise(self):
    
        """
        
        Method that finishes building the piecewise eos by computing the three core tropic coefficients and instantiating a 
        Polytropic object for each layer. It also computes transition energy densitites and pressures
        -------

        """
    
        self.build_k()
        
        for k, gamma in zip(self.kappas, self.gammas):           
            self.layers.append(Polytropic(k, gamma))
            
        self.densities[4] = (self.kappas[4]/self.kappas[3])**(1/(self.gammas[3]-self.gammas[4]))
        
        prev_layer = None
        for (layer, density) in zip(self.layers, self.densities):

                if not(prev_layer == None):
                    layer.a = self.set_transition(prev_layer, layer, density)
                else:
                    density = 0.0

                eden = layer.eden_from_density(density)
                pressure = layer.pressure_from_density(density)

                self.pressures.append(pressure)
                self.edens.append(eden)
                
                prev_layer = layer



    def set_transition(self, prev_layer, layer, density):
    
        """
        
        Method that computes the transition continuity constant of a layer given the current layer and the previous one

        Parameters
        ----------
        prev_layer : Polytropic object
            previous layer
        layer : Polytropic object
            current layer
        density : float
            density of transition between layers

        Returns transition : float
            a + (k/(gamma-1))*d^(gamma-1) of the previous layer - (k/(gamma-1))*d^(gamma-1) of the current layer
        -------

        """   
        
        transition = prev_layer.a + (prev_layer.k/(prev_layer.gamma - 1))*density**(prev_layer.gamma-1) - (layer.k/(layer.gamma - 1))*density**(layer.gamma-1)
        return transition
        
        
        
    def pressure_from_density(self, density):  
        
        """
        
        Method that returns the pressure given the density, by finding the layer the density belongs to, and then returning 
        the polytropic pressure

        Parameters
        ----------
        density : float
            input density

        Returns pressure : float
            computed pressure
        -------

        """
        
        layer = self.find_layer(density, "density")
        pressure = layer.pressure_from_density(density) 
        return pressure



    def density_from_pressure(self, pressure):  

        """
        
        Method that returns the density given the pressure, by finding the layer the pressure belongs to, and then returning 
        the polytropic density

        Parameters
        ----------
        pressure : float
            input pressure

        Returns density : float
            computed density
        -------

        """ 
             
        layer = self.find_layer(pressure, "pressure")
        density = layer.density_from_pressure(pressure) 
        return density



    def eden_from_pressure(self, pressure): 

        """
        
        Method that returns the energy density given the pressure, by finding the layer the pressure belongs to, and then returning 
        the polytropic energy density

        Parameters
        ----------
        pressure : float
            input pressure

        Returns eden : float
            computed energy density
        -------

        """
               
        layer = self.find_layer(pressure, "pressure")
        eden = layer.eden_from_pressure(pressure)   
        return eden
    
    
    
    def find_layer(self, value, value_type):
    
        """
        
        Method that finds a layer given a value of pressure or density, by confronting it with the transition
        densities or pressures of the Piecewise

        Parameters
        ----------
        value : float
            input pressure or density
        value_type : string
            "pressure" or "density"

        Returns layers[i] : Polytropic object
            the found layer 
        -------

        """
        
        if value_type == "pressure":
            
            if value >= self.pressures[-1]:
                return self.layers[-1]
            
            for i in range( len(self.pressures) - 1):
                if self.pressures[i] <= value < self.pressures[i+1]:
                    return self.layers[i]
                
        elif value_type == "density":
            
            if value >= self.densities[-1]:
                return self.layers[-1]
            
            for i in range( len(self.densities) - 1):
                if self.densities[i] <= value < self.densities[i+1]:
                    return self.layers[i] 
                
                
                
class Implicit:
    
    
    
    def __init__(self):
    
        """
        
        Class that, given the explicit expression of the energy and the pressure in function of x = k_f/(m*c), where k_f 
        is the Fermi momentum, returns the interpolated eos e = e(p)
        -------

        """
        
        self.kind = "PressureEnergyDensityImplicit"        
        self.e0 = ((ut.M_N**4)*(ut.C_SI**5)/((np.pi**2)*(ut.HBAR_SI**3)))*cgs_geom_dictionary["si"]["energy_density"]["geom"]
        self.density_from_pressure, self.pressure_from_density = self.interpolate_solution()
        self.eden_from_pressure = self.density_from_pressure
        
        
        
    def implicit_density(self, x):
    
        """
        
        Method that returns the density corresponding to a reduced k_fermi x
        
        Parameters
        ----------
        x : float
            k_fermi/m*c
        -------

        """
        
        imp_density = (self.e0/8)*((2*x**3 + x)*np.sqrt(1.0 + x**2) - np.arcsinh(x))
        return imp_density
    


    def implicit_pressure(self, x):
        
        """
        
        Method that returns the pressure corresponding to a reduced k_fermi x
        
        Parameters
        ----------
        x : float
            k_fermi/m*c
        -------

        """
        
        imp_pressure = (self.e0/24)*((2*x**3 - 3*x)*np.sqrt(1.0 + x**2) + 3*np.arcsinh(x))
        return imp_pressure
    
    
    
    def interpolate_solution(self):
    
        """
        
        Method that interpolates two array of 100000 pressures and energy densities

        Returns cubic_spline : function
            return of utils.cubic_spline
        -------

        """
        
        x_range = np.linspace(0, 1e3, int(100000))
        p_column, e_column = self.implicit_pressure(x_range), self.implicit_density(x_range)

        return cubic_spline(p_column, e_column), cubic_spline(e_column, p_column)

        

    
    
    
    
    
    
    
    
    
        
        
        
        