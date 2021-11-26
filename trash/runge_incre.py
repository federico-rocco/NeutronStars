# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 16:33:15 2021

@author: cosmo
"""

def step_comp(f, u, t, dt):
    
    g0 = f(u, t)    
    g1 = f(u + 0.5*dt*g0, t+0.5*dt)    
    g2 = f(u + 0.5*dt*g1, t+0.5*dt)    
    g3 = f(u + dt*g2, t+dt)
        
    return (1/6)*dt*(g0 + 2*g1 + 2*g2 + g3)
    
def rk4(f, u, t, dt, csi):
        
    y1 = step_comp(f, u, t, dt)   
    y2 = step_comp(f, u, t, dt/2)
    err = abs((y1 - y2)/15) 
    dt = dt*(csi/err)**(0.2)
        
    return y1