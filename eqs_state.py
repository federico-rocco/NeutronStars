# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 12:13:14 2021

@author: cosmo
"""

class eqs_state:
    
    eos = ""
    rs = 0
    
    #def __init__(self, kind, rho_s):
        #if kind != "a","b","c": mettere error
        #self.eos = kind
        #self.rs = rho_s
      
    def select():
        eqs_state.eos = input("Che eq voui? a: b: c: ") 
        eqs_state.rs = float(input("Che valore iniziale voui?"))
        #print(eqs_state.rs)
        
    def density_calc(pres):
        if eqs_state.eos == "prima":
            n = (pres/363.44)**(1/2.54)
            rho = 236*n**2.54+n*938
        
        elif eqs_state.eos == "seconda":
            n = 4*pres
            rho = 45*n
            
        return rho
    
    def central_density():
        #print(eqs_state.rs)
        return eqs_state.rs
        