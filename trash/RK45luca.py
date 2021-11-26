# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 11:39:34 2021

@author: lucab
"""

#RK45 prova
#gli array non vanno bene per ora perché hanno dimensione fixata e non posso
#aggiungere argomenti a ruota. Fare questo mi serve per avere un grafico decente
#e due array decenti. Non so se mi serve però
import numpy as np 
import matplotlib.pyplot as plt 

a = 0
b= 2*np.pi
L = b - a
csi = 10**(-10)
N = 200

U_0 = 0
t_0 = a

t = np.zeros(N)
u = np.zeros(N)
z =  np.zeros(N)

t[0] = t_0
u[0] = U_0
z[0] = U_0


def f(x, y):
    return np.cos(x)

for i in range(N-1): 
    d = 0
    h = 1
    e = 0.01
    while d < 1:
        k1 = h*f(t[i], u[i])
        k2 = h*f(t[i] + h/4, u[i] + k1/4)
        k3 = h*f(t[i] + 3*h/8, u[i] + 3*k1/32 + 9*k2/32)
        k4 = h*f(t[i] + 12*h/13, u[i] + 1932*k1/2197 - 7200*k2/2197 + 7296*k3/2197)
        k5 = h*f(t[i] + h, u[i] + 439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104)
        k6 = h*f(t[i] + h/2, u[i] - 8*k1/27 + 2*k2 -3544*k3/2565 + 1859*k4/4104 -11*k5/40)
        u[i+1] = u[i] + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5
        z[i+1] = u[i] + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 + 2*k6/55
        
        R = np.abs(z[i+1] - u[i+1])
        d = 0.9*(csi/R)**(0.25)
        h = h*d
    print(d)
    e = e*d
    k1 = e*f(t[i], u[i])
    k2 = e*f(t[i] + e/4, u[i] + k1/4)
    k3 = e*f(t[i] + 3*e/8, u[i] + 3*k1/32 + 9*k2/32)
    k4 = e*f(t[i] + 12*e/13, u[i] + 1932*k1/2197 - 7200*k2/2197 + 7296*k3/2197)
    k5 = e*f(t[i] + e, u[i] + 439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104)
    k6 = e*f(t[i] + e/2, u[i] - 8*k1/27 + 2*k2 -3544*k3/2565 + 1859*k4/4104 -11*k5/40)
    u[i+1] = u[i] + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5
    t[i+1] = t[i] + e

plt.plot(t, u)
plt.show()