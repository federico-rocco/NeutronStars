# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 16:09:43 2021

@author: cosmo
"""
import numpy as np
from scipy.constants import speed_of_light, gravitational_constant, hbar, m_n
import astropy.constants


C_SI = speed_of_light #m/s
G_SI = gravitational_constant  #m^3/(kg*s^2)
MSUN_SI = astropy.constants.M_sun.value  #kg
MSUN_CGS = MSUN_SI*1000 #g
HBAR_SI = hbar

cgs_geom_dictionary = { "geom": { "lenght": {"cm": 100.,
                                             "m": 1,
                                             "km": 1e-3
                                             },
                                 
                                  "time": {"cgs": 1/C_SI,
                                           "si": 1/C_SI,
                                           "geom": 1.
                                           },
                                 
                                  "mass": { "g": 1e3*(C_SI ** 2.) / G_SI,
                                            "kg": (C_SI ** 2.) / G_SI,
                                            "m_sol": (C_SI ** 2.) / G_SI / MSUN_SI
                                           },
                                 
                                  "density": {"cgs": 1e-3*(C_SI ** 2.) / G_SI
                                              },
                                  
                                  "energy": {"cgs": 1e7*(C_SI ** 4.) / G_SI, 
                                             "si": (C_SI ** 4.) / G_SI,
                                             },
                                  
                                  "energy_density": {"cgs": 10*(C_SI ** 4.) / G_SI
                                                     },
                                  
                                  "pressure": {"cgs": 10*(C_SI ** 4.) / G_SI
                                               }
                                  },
                       
                        "cgs": { "lenght": { "cm": 1.,
                                            "m": 1e-2,
                                            "km": 1e-5
                                            },
                                
                                 "time": { "cgs": 1,
                                          "si": 1,
                                          "geom": C_SI
                                          },
                                
                                 "mass": {"kg": 1e-3,
                                          "geom": 1e-3*G_SI/(C_SI ** 2.) ,
                                          "m_sol": 1/MSUN_CGS
                                          },
                                                        
                                 "density": {"si": 1e3,
                                             "geom": 1e3*G_SI/(C_SI ** 2.)
                                             },
                                                         
                                 "energy": {"si": 1e-7, 
                                            "geom": 1e-7*G_SI/(C_SI ** 4.)
                                            },
                                                         
                                 "energy_density": {"si": 0.1,
                                                    "geom": 0.1*G_SI/(C_SI ** 4.) 
                                                    },
                                                         
                                 "pressure": {"si": 0.1 ,
                                              "geom": 0.1*G_SI/(C_SI ** 4.)
                                              }
                                 },
                        
                        "si" : { "lenght": {"cm": 100.,
                                            "m": 1.,
                                            "km": 1e-3
                                            },
                                
                                 "time": {"cgs": 1.,
                                          "geom": C_SI,
                                          },
                                 
                                 "mass": {"g": 1000.,
                                          "kg": 1.,
                                          "geom":  G_SI/(C_SI ** 2.),
                                          "m_sol": 1 / MSUN_SI
                                          },
                                
                                 "density": {"cgs": 1e-3,
                                             "geom": (C_SI ** 2.) / G_SI,
                                             },
                                                                
                                 "energy": {"cgs": 1e7, 
                                            "geom":  G_SI/(C_SI ** 4.)
                                           },
                                
                                "energy_density": {"cgs": 10.,
                                                   "geom": G_SI/(C_SI ** 4.)
                                                   },
                                
                                "pressure": {"cgs": 10.,
                                             "geom": G_SI/(C_SI ** 4.)
                                             },
                                }
                       }


def step_comp(f, t, dt, u, v):
        
    k1 = f(t, [u, v])
    k2 = f(t + dt/4, [u + dt*k1[0]/4, v + dt*k1[1]/4])
    #print("k1",k1)
    #print("v=",v,"k11=",k1[1],"k21=",k2[1],"dt=",dt,"v + 3*dt*k1[1]/32 + 9*dt*k2[1]/32]=",v + 3*dt*k1[1]/32 + 9*dt*k2[1]/32)
    k3 = f(t + 3*dt/8, [u + 3*dt*k1[0]/32 + 9*dt*k2[0]/32, v + 3*dt*k1[1]/32 + 9*dt*k2[1]/32])
    k4 = f(t + 12*dt/13, [u + 1932*dt*k1[0]/2197 - 7200*dt*k2[0]/2197 + 7296*dt*k3[0]/2197, v + 1932*dt*k1[1]/2197 - 7200*dt*k2[1]/2197 + 7296*dt*k3[1]/2197])
    k5 = f(t + dt, [u + 439*dt*k1[0]/216 - 8*dt*k2[0] + 3680*dt*k3[0]/513 - 845*dt*k4[0]/4104, v + 439*dt*k1[1]/216 - 8*dt*k2[1] + 3680*dt*k3[1]/513 - 845*dt*k4[1]/4104])
    k6 = f(t + dt/2, [u - 8*dt*k1[0]/27 + 2*dt*k2[0] -3544*dt*k3[0]/2565 + 1859*dt*k4[0]/4104 -11*dt*k5[0]/40, v - 8*dt*k1[1]/27 + 2*dt*k2[1] -3544*dt*k3[1]/2565 + 1859*dt*k4[1]/4104 -11*dt*k5[1]/40])
    #print("k2=",k2,"k3=",k3) 
    #print(k1,k2,k3,k4,k5)
    du = 25*dt*k1[0]/216 + 1408*dt*k3[0]/2565 + 2197*dt*k4[0]/4101 - dt*k5[0]/5
    dv = 25*dt*k1[1]/216 + 1408*dt*k3[1]/2565 + 2197*dt*k4[1]/4101 - dt*k5[1]/5
    #print("dm=",du,"dp=",dv)
    du1 = 16*dt*k1[0]/135 + 6656*dt*k3[0]/12825 + 28561*dt*k4[0]/56430 - 9*dt*k5[0]/50 + 2*dt*k6[0]/55
    dv1 = 16*dt*k1[1]/135 + 6656*dt*k3[1]/12825 + 28561*dt*k4[1]/56430 - 9*dt*k5[1]/50 + 2*dt*k6[1]/55
        
    return [du, dv,du1,dv1]
    
def adaptiveRungeKutta(f, t, u, v, dt, csi, hmax):
    
    y = step_comp(f, t, dt, u, v)

    err = max(abs(y[1] - y[3]), abs(y[1] - y[3]))
    if err != 0:
        new_step = 0.9*dt*(csi*dt/err)**(0.25)
        if new_step < csi:
            dt = csi
        elif new_step > hmax:
            dt = hmax
        else:
            if new_step < dt:
                dt = new_step
                y = step_comp(f, t, dt, u, v) 
            else:
                dt = new_step
        
    return [y[0], y[1], dt]


def ODEsolver(eq_type, central_pressure, eq_state):
      
    step = 0.1
    csi = 1e-5
    max_step = 10
    
    mass = np.array([0])
    pressure = np.array([central_pressure])
    radius = np.array([csi])
    
    i=0
    while(pressure[i]>0):
        
        dy = adaptiveRungeKutta(eq_type, radius[i], mass[i], pressure[i], step, csi, max_step)
        step = dy[2]
        
        if not pressure[i]+dy[1] > 0:
            break
        mass = np.append(mass, mass[i] + dy[0])
        pressure = np.append(pressure, pressure[i] + dy[1])
        radius = np.append(radius, radius[i]+step)          
        
        i = i+1

        #print("step:",step,"radius:",radius[i],"mass:",mass[i],"pressure:",pressure[i])

    return radius,mass,pressure


import bisect

def CubicSpline(x, y):
#costruisci gli intervalli dal tuo array di dati
    n = len(x)
    h = []
    for i in range(n-1):
        h.append(x[i+1] - x[i])
#per risolvere il sistema di equazioni della spline conviene usare la derivata
#seconda, che per una cubica è lineare in x. Integrando due volte e usando la
#continuità togliamo le 2n costanti di integrazione. Poi derivando e usando
#il matching della derivata prima esplicitiamo n-1 M. I restanti vengono dalle
#boundary che sono settate a zero (Boundary type 2).
#Il sistema da risolvere risulta fatto da una matrice triangolare e per risolverlo
#introduciamo c_p e d_p. 
    A = []
    B = [] 
    C = [] 
    C.append(0)
    for i in range(n-2):
        A.append(h[i]/(h[i] + h[i+1])) 
        C.append(h[i+1]/(h[i] + h[i+1]))        
    for i in range(n):
       B.append(2)
    A.append(0)
#RHS del sistema di equazioni con la tridiagonal
    D = [] 
    D.append(0)
    for i in range(1, n-1):
        D.append(6*((y[i+1] - y[i])/h[i] - (y[i] - y[i-1])/h[i-1])/(h[i] + h[i-1]))
    D.append(0)
#soluzioni del sistema con metodo di Thomas 
    c_p = C + [0] 
    d_p = [0]*n
    M = [0]*n 
    c_p[0] = C[0]/B[0]
    d_p[0] = D[0]/B[0]
    for i in range(1, n):
        c_p[i] = (c_p[i]/(B[i] - c_p[i-1]*A[i-1]))
        d_p[i] = (D[i] - d_p[i-1]*A[i-1])/(B[i] - c_p[i-1]*A[i-1])
    """clean"""
#M è la derivata seconda della spline in tutti i punti, soluzione del sistema sopra
    M[-1] = d_p[-1]
    for i in range(n-2, -1, -1):
        M[i] = d_p[i] - c_p[i]*M[i + 1]
    
    #formula per i coefficienti in una cubic spline 
    coefficients = []
    for i in range(n-1):
        coefficients.append([((M[i+1] - M[i])*h[i]**2)/6, (M[i]*h[i]**2)/2, (y[i+1] - y[i] - ((M[i+1] + 2*M[i])*h[i]**2)/6), y[i]])
#bisect(list, num):This function returns the position in the sorted list, 
#where the number passed in argument can be placed so as to maintain the resultant list in sorted order. 
#If the element is already present in the list, the right most position where element has to be inserted is returned.
    def spline(val):    
        idx = min(bisect.bisect(x, val) - 1, n - 2)
        z = (val - x[idx])/h[idx]
        Coef = coefficients[idx]
        S_z = Coef[0]*z**3 + Coef[1]*z**2 + Coef[2]*z + Coef[3]
        return S_z 
    return spline 