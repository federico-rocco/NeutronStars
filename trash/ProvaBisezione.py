# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 18:48:45 2021

@author: lucab
"""


import numpy as np
import matplotlib 
import matplotlib.pyplot as plt    
    
def pres(x):
    mn =  938    
    e0 = mn**4/(np.pi**2)
    return (e0/24)*((2*x**3 - 3*x)*(1 + x**2)**0.5 + 3*np.arcsinh(x))


def bisection(func, p, interval):

    tol = 0.01    
    a, b = interval[0], interval[1]           
    kfermi = np.array([])
    
    j = 0
    for i in range(p.size):
        f_a = func(a)-p[i]
        while b - a > tol :
            m = (a + b)/2
            f_m = func(m) - p[i]  
               
            if np.sign(f_a) == np.sign(f_m):
                a= m
                f_a = f_m
            else:
                b = m
            j += j
            
        kfermi = np.append(kfermi, m) 

    return kfermi

p = np.linspace(10**10, 10**50, 20)
kfermi_range = np.array([-1000000000000000000000000000000000000, 100000000000000000000000000000000000000000000000000000000000000])
kfermi_range = kfermi_range.astype(float)
k = bisection(pres, p, kfermi_range)

plt.plot(k, p)
plt.show()


    
