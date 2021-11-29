# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 16:09:43 2021

@author: cosmo
"""
import numpy as np



def step_comp(f, t, dt, u, v):
        
    k1 = f(t, u, v)
    k2 = f(t + dt/4, u + dt*k1[0]/4, v + dt*k1[1]/4)
    k3 = f(t + 3*dt/8, u + 3*dt*k1[0]/32 + 9*dt*k2[0]/32, v + 3*dt*k1[1]/32 + 9*dt*k2[1]/32)
    k4 = f(t + 12*dt/13, u + 1932*dt*k1[0]/2197 - 7200*dt*k2[0]/2197 + 7296*dt*k3[0]/2197, v + 1932*dt*k1[1]/2197 - 7200*dt*k2[1]/2197 + 7296*dt*k3[1]/2197)
    k5 = f(t + dt, u + 439*dt*k1[0]/216 - 8*dt*k2[0] + 3680*dt*k3[0]/513 - 845*dt*k4[0]/4104, v + 439*dt*k1[1]/216 - 8*dt*k2[1] + 3680*dt*k3[1]/513 - 845*dt*k4[1]/4104)
    k6 = f(t + dt/2, u - 8*dt*k1[0]/27 + 2*dt*k2[0] -3544*dt*k3[0]/2565 + 1859*dt*k4[0]/4104 -11*dt*k5[0]/40, v - 8*dt*k1[1]/27 + 2*dt*k2[1] -3544*dt*k3[1]/2565 + 1859*dt*k4[1]/4104 -11*dt*k5[1]/40)
        
    du = 25*k1[0]/216 + 1408*k3[0]/2565 + 2197*k4[0]/4104 - k5[0]/5
    dv = 25*k1[1]/216 + 1408*k3[1]/2565 + 2197*k4[1]/4104 - k5[1]/5
    return [du, dv]
    
def rk4(f, t, u, v, dt, csi):
    
    y1 = step_comp(f, t, dt, u, v) 
    y2 = step_comp(f, t, dt/2, u, v)
    err = max(abs((y1[0] - y2[0])/15), abs((y1[0] - y2[0])/15)) 
    if err != 0:
        new_step = dt*(csi/err)**(0.2)
        if new_step < dt:
            dt = new_step
            y1 = step_comp(f, t, dt, u, v) 
                
        else:
            dt = new_step
        
    return y1

def solve(eq_type, central_value, eq_state):
      
    rho = 1
    #print(rho)
    step = 0.001
    csi = 0.001
    G = 6.67*10**(-11)
    
    mass = np.array([0])
    pres = np.array([rho])
    radius = np.array([csi])
    
    """if eq_type == "Newton":
        func = ns.Newton_eqs
    elif eq_type == "TOV":
        func = ns.TOV_eqs"""
    
    i=0
    while(pres[i] > 0):
        print("Ciclo ", i)
        dy = rk4(eq_type, radius[i], mass[i], pres[i], step, csi)
        mass = np.append(mass, mass[i] + dy[0])
        pres = np.append(pres, pres[i] + dy[1])
        radius = np.append(radius, radius[i]+step)
        rho = eq_state.energy_density(pres[i])
        
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

def bisection(f, xL, xR, guess):
    fL = f(xL)- guess
    xM = (xL + xR)/2.0
    fM = f(xM) - guess
    epsilon = 10**(-10)
    while abs(fM) > epsilon:
        xM = (xL + xR)/2.0
        if fL*fM > 0: 
            xL = xM
            fL = fM
        else:
            xR = xM
        xM = (xL + xR)/2.0
        fM = f(xM) - guess
    return xM
