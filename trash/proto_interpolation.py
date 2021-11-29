# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 17:21:30 2021

@author: lucab
"""

import numpy as np 
from scipy import interpolate
import matplotlib.pyplot as plt 

def energy(x):
    return (1/8)*((2*x**3 + x)*np.sqrt(1 + x**2) - np.arcsinh(x))

def pressure(x):
    return (1/24)*((2*x**3 - 3*x)*np.sqrt(1 + x**2) + 3*np.arcsinh(x))
## va sistemato il prefactor e0/8##

a=-250
b=250
mn =  1.6*10**(-30) 
h_bar = 6*10**(-27)
c = 3*10**10
e0 = ((mn**4)*(c**5))/(((np.pi)**2)*(h_bar)**3)
e_guess = np.linspace(0, b-a, b-a)
eps = 10**(-10)

def bisection(f, xL, xR, guess):
    fL = f(xL)- guess
    xM = float((xL + xR)/2.0)
    fM = f(xM) - guess
    while abs(fM) > eps:
        xM = (xL + xR)/2.0
        if fL*fM > 0: 
            xL = xM
            fL = fM
        else:
            xR = xM
        xM = float((xL + xR)/2.0)
        fM = f(xM) - guess
    return xM

sol = []

for i in range(b-a):
    sol.append(bisection(energy, a, b, e_guess[i]))

eos = []

for j in range(b-a):
    eos.append(pressure(sol[j]))

rel_eos = interpolate.interp1d(e_guess, eos, kind = 3)
print(rel_eos(0.5))

"""p = np.zeros(b-a)

for k in range(b-a):
    p[k]=(e_guess[k])**(4/3)    
plt.plot(e_guess, eos, e_guess, p)

plt.show()"""


        
