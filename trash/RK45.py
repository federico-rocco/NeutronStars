# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 11:12:45 2021

@author: lucab
"""

#RK45 

import numpy as np

import matplotlib.pyplot as plt

a = 0

b = 10

L = b - a

h = 0.1

csi = 10**(-10)

N = int(L/h)

u = np.zeros(N)

t = np.linspace(a, b, N)

u0 = 1

def RK45(f, U_0, alfa, tol):
    
    z = np.zeros(N)
    
    z[0] = U_0
    
    u[0] = U_0
    
    for n in range(N-1):
        
        k1 = f(t[n], u[n])
        
        k2 = f(t[n] + alfa/4, u[n] + k1/4)
        
        k3 = f(t[n] + 3*alfa/8, u[n] + 3*k1/32 + 9*k2/32)
        
        k4 = f(t[n] + 12*alfa/13, u[n] + 1932*k1/2197 + 7200*k2/2197 + 7296*k3/2197)
        
        k5 = f(t[n] + alfa, u[n] + 439*k1/216 - 8*k2 + 3680*k3/513 + 845*k4/4104)
        
        k6 = f(t[n] + 0.5*alfa, u[n] -8*k1/27 + 2*k2 - 3544*k2/2565 + 1859*k4/4104 - 11*k5/40)
        
        u[n+1] = u[n] + alfa*(25*k1/216 + 1408*k3/2565 + 2197*k4/4101 - k5/5)
        
        z[n+1] = u[n] + alfa*(16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55)
        
        beta = (tol*alfa/(2*abs(z[n+1] - u[n+1])))**(0.25)
        
        alfa = alfa*beta
        
    return u 

def demo(x, y):
    
    return -2*y + np.exp(-2*(x - 6)**2)

u = RK45(demo, u0, h, csi)

plt.plot(t, u)

plt.show()
        
        
        
        
        
        
        
        
        
        
        
        