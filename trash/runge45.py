# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 14:32:59 2021

@author: cosmo
"""

#RK45 

import numpy as np
import matplotlib.pyplot as plt

def RK45(f, U_0, step, csi):  
    z = np.zeros(N)
    z[0] = U_0  
    u[0] = U_0
    for n in range(N-1):
        print("step: ", step)
        k1 = step*f(t[n], u[n])
        k2 = step*f(t[n] + step/4, u[n] + k1/4)
        k3 = step*f(t[n] + 3*step/8, u[n] + 3*k1/32 + 9*k2/32)
        k4 = step*f(t[n] + 12*step/13, u[n] + 1932*k1/2197 - 7200*k2/2197 + 7296*k3/2197)
        k5 = step*f(t[n] + step, u[n] + 439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104)
        k6 = step*f(t[n] + 0.5*step, u[n] -8*k1/27 + 2*k2 - 3544*k3/2565 + 1859*k4/4104 - 11*k5/40) 
        u[n+1] = u[n] + 47*k1/450 + 12*k3/25 + 32*k4/225 + k5/30 +6*k6/25
        #z[n+1] = u[n] + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55
        print("u: ", u[n+1], " z: ", z[n+1])
        err = abs(-k1/150+3*k3/100-16*k4/75-k5/20+6*k6/25)
        print(err)
        if err != 0:
            step = 0.9*step*(csi/(err))**(0.25)   
    return u 

def demo(x, y): 
    return np.cos(x)

a = 0
b = 10
h = 0.1
csi = 10**(-5)
N = int((b-a)/h)
u = np.zeros(N)
t = np.linspace(a, b, N)
u0 = 0

u = RK45(demo, u0, h, csi)
plt.plot(t, u)
plt.show()
        
        
        
        
        
        
        
        
        
        
        
        