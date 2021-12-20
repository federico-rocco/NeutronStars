# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 19:02:16 2021

@author: cosmo
"""

import numpy as np
import multiprocessing as mp
from utils import ode_solver, cgs_geom_dictionary

class NeutronStar:

    
    def __init__(self, NS_type, eq_state):
    
        """
        
        Class that, given a specific equation of state, returns a neutron star with the appropriate structure

        Parameters
        ----------
        NS_type : string
            type of neutron star, e.g. Non Relativistic Pure NS, Relativistic Pure NS or Non Pure NS
            it keeps track of the kind of star in the case different stars are needed
        eq_state : object of eqs_state module (Polytropic, Piecewise or Implicit)
            the equation of state that the star needs to follow, linking the energy density and the pressure in each point of the star
            Non Relativistic Pure NS -> Polytropic
            Relativistic Pure NS -> Implicit
            Non Pure NS -> Piecewise
        -------

        """
        
        self.kind = NS_type
        self.eos = eq_state
        
    
    def tov_eqs(self, r, y):
        
        """
        
        Tolman-Oppenheimer-Volkoff equations: differential equations that describe the variation of mass and pressure inside the star
        They add relativistic corrections to the Newton equations

        Parameters
        ----------
        r : float
            independent variable, it represents the distance from the centre of the star
        y : array of float
            dependent variable, array containing the mass m = y[0] and the pressure p = y[1]

        Returns dy : array of float containing the increments of mass and pressure
        -------
        Note: equations written in the form dy/dx = f(x,y)
        eden is the energy density corresponding to the pressure p, computed using EdenFromPressure of the eos
        """
        
        m, p = y
        if p < 0:
            return [0,0]
        eden = self.eos.eden_from_pressure(p)
        dm = 4*np.pi*eden*r**2
        dp = - (eden+p)*(m + 4*np.pi*r**3*p)/(r*(r - 2*m))
        dy = [dm, dp]
        
        return dy

    def newton_eqs(self, r, y):

        """
        
        Newton equations: differential equations that describe the variation of mass and pressure inside the star
        They are hydrostatic equilibrium equations and provide the non-relativistic limit of the TOV equations

        Parameters
        ----------
        r : float
            independent variable, it represents the distance from the centre of the star
        y : array of float
            dependent variable, array containing the mass m = y[0] and the pressure p = y[1]

        Returns dy : array of float containing the increments of mass and pressure
        -------
        Note: equations written in the form dy/dx = f(x,y)
        eden is the energy density corresponding to the pressure p, computed using EdenFromPressure of the eos
        """
        m, p = y
        if p < 0:
            return [0,0]
        eden = self.eos.eden_from_pressure(p)
        dm = 4*np.pi*eden*r**2
        dp = - (eden*m)/(r**2)
        dy = [dm, dp]

        return dy
    
    def star_solver(self, eq_type, central_value, value_type="pressure", unit_type="geom"):

        """
        
        Method that solves the structure of the star by solving a system of differential equations        

        Parameters
        ----------
        eq_type : newton_eqs or tov_eqs, methods of NeutronStar class
            allows to choose bewtween Newton equations and TOV equations
        central_value : float
            value of the pressure or density in the centre of the star
            it represents one of the Cauchy conditions to solve the differential system
        value_type : string, optional
            it allows to pass a pressure or a density as central value. The default is "pressure". If "density" is chosen the value
            is converted into a pressure
        unit_type : string, optional
            it specifies the system of measurement in which central_value is expressed. The default is "geom". If "cgs" or "si" are 
            chosen, the value is converted accordingly

        Returns r_out, m_out, p_out, arrays of float containing the radii, masses and pressures computed in each iteration. The last
        elements of each array correspond to the total radius, the total mass and the last pressure (less than zero)
        -------
        Note: the solver is set to work in geometric units in order to reduce numerical errors

        """
        

        if eq_type not in [self.tov_eqs, self.newton_eqs]:
            raise ValueError("Differential equations to be passed: tov_eqs os newton_eqs") 
        if value_type not in ["pressure", "density"]:
            raise ValueError("Can only pass a pressure or a density as central value") 
        if unit_type not in ["geom", "si", "cgs"]:
            raise ValueError("Specify \"geom\", \"si\" or \"cgs\" as unit type") 
  
        central_value = central_value*cgs_geom_dictionary[unit_type][value_type]["geom"]
   
        if value_type == "density":
            central_value = self.eos.pressure_from_density(central_value)

        solutions = ode_solver(eq_type, central_value)
        r_out = np.array([])
        m_out = np.array([])
        p_out = np.array([])

        for i in range(solutions[0].size):
            r_out = np.append(r_out, solutions[0][i]*cgs_geom_dictionary["geom"]["lenght"]["km"])
            m_out = np.append(m_out, solutions[1][i]*cgs_geom_dictionary["geom"]["mass"]["m_sol"])
            p_out = np.append(p_out, solutions[2][i]*cgs_geom_dictionary["geom"]["pressure"]["cgs"])

        return r_out, m_out, p_out
        
    
    def mass_vs_radius_old(self, centrals, eq_type, value_type, unit_type):
        
        radii = np.array([])
        masses = np.array([])

        
        import tqdm               
        for i in tqdm.tqdm(range(centrals.size)):
            solutions = self.star_solver(eq_type, centrals[i], value_type, unit_type)
            radii = np.append(radii, solutions[0][-1])
            masses = np.append(masses, solutions[1][-1])
        print(__name__)    
        return radii, masses
    
    
    def mass_vs_radius(self, centrals, eq_type, value_type="pressure", unit_type="geom"):
        
        """
        
        Method that solves a sequence of stars of the same type by calling star_solver. Multiprocessing is used to reduce 
        computational time
        
        Parameters
        ----------
        centrals : array of float
            array containing the central value of each star
        eq_type : Newton_eqs or TOV_eqs, methods of NeutronStar class
            allows to choose bewtween Newton equations and TOV equations
        value_type : string, optional
            it allows to pass a pressure or a density as central value. The default is "pressure"
        unit_type : string, optional
            it specifies the system of measurement in which central_value is expressed. The default is "geom"

        Returns radii and masses, arrays of float containing the total radius and mass of all the stars in the sequence
        -------

        """ 
        
        radii = np.array([])
        masses = np.array([])
            
        import tqdm
        pool = mp.Pool(mp.cpu_count())
        for _ in tqdm.tqdm(pool.apply(self.star_solver, args=(eq_type, central, value_type, unit_type)), total=len(centrals)):
            pass
        results = [pool.apply(self.star_solver, args=(eq_type, central, value_type, unit_type)) for central in centrals]
        pool.close()
        pool.join()
    
        for result in results:
            radii = np.append(radii, result[0][-1])
            masses = np.append(masses, result[1][-1])
            
        return radii, masses
        
            
            