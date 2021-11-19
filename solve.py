# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 16:09:43 2021

@author: cosmo
"""
def solve():
    import numpy as np
    from project.eqs_state import eqs_state
    
    def step_comp(f, u, t, dt, rho, m=""):
        
        g0 = f(u, t, rho, m)    
        g1 = f(u + 0.5*dt*g0, t+0.5*dt, rho, m)    
        g2 = f(u + 0.5*dt*g1, t+0.5*dt, rho, m)    
        g3 = f(u + dt*g2, t+dt, rho, m)
            
        return (1/6)*dt*(g0 + 2*g1 + 2*g2 + g3)
        
    def rk4(f, u, t, dt, csi, rho, m=""):
        
        y1 = step_comp(f, u, t, dt, rho, m) 
        if f==func2:
            y2 = step_comp(f, u, t, dt/2, rho, m)
            err = abs((y1 - y2)/15) 
            if err != 0:
                dt = dt*(csi/err)**(0.2)
            
        return y1
    
    
    rho = eqs_state.central_density()
    #print(rho)
    step = 0.001
    csi = 0.001
    G = 6.67*10**(-11)
    
    mass = np.array([0])
    pres = np.array([rho])
    radius = np.array([csi])
    
    
    i=0
    while(pres[i] > 0):
        print("Ciclo ", i)
        mass = np.append(mass, mass[i] + rk4(func1, mass[i], radius[i], step, csi, rho))
        pres = np.append(pres, pres[i] + rk4(func2, pres[i], radius[i], step, csi, rho, mass[i]))
        radius = np.append(radius, radius[i]+step)
        rho = eqs_state.density_calc(pres[i])
        
        if pres[i+1] < 0:
            #print("R= ", radius[i])
            break
        else:
            print("radius:", radius[i+1],"mass:", mass[i+1], "pressure=",pres[i+1])
        i+=1
    
    #print(rho)
    mass = mass*4*np.pi*rho/(4*np.pi*rho*G)**1.5
    radius = radius/(np.sqrt(4*np.pi*rho*G))
    print("Total mass: ", mass[mass.size-1], "Radius: ", radius[radius.size-1])
    
    return radius,mass,pres

def func1(y, r, rho, mass):
    return rho*r**2
    
def func2(y, r, rho, mass):
    return -mass*rho*r**(-2)