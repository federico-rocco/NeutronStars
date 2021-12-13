# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 16:09:43 2021

@author: cosmo
"""
import numpy as np
from scipy.constants import speed_of_light, gravitational_constant, hbar, m_n
import astropy.constants
import bisect


C_SI = speed_of_light #m/s
C_CGS = C_SI*100
G_SI = gravitational_constant  #m^3/(kg*s^2)
G_CGS = G_SI*1000
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
        
    k1 = f(t, [u, v],True)
    
    if v + dt*k1[1]/4<0:
        return [0,0,0,0,False]
    k2 = f(t + dt/4, [u + dt*k1[0]/4, v + dt*k1[1]/4])
    
    if v + 3*dt*k1[1]/32 + 9*dt*k2[1]/32<0:
        return [0,0,0,0,False]
    k3 = f(t + 3*dt/8, [u + 3*dt*k1[0]/32 + 9*dt*k2[0]/32, v + 3*dt*k1[1]/32 + 9*dt*k2[1]/32])
    
    if v + 1932*dt*k1[1]/2197 - 7200*dt*k2[1]/2197 + 7296*dt*k3[1]/2197<0:
        return [0,0,0,0,False]
    k4 = f(t + 12*dt/13, [u + 1932*dt*k1[0]/2197 - 7200*dt*k2[0]/2197 + 7296*dt*k3[0]/2197, v + 1932*dt*k1[1]/2197 - 7200*dt*k2[1]/2197 + 7296*dt*k3[1]/2197])
    
    if v + 439*dt*k1[1]/216 - 8*dt*k2[1] + 3680*dt*k3[1]/513 - 845*dt*k4[1]/4104<0:
        return [0,0,0,0,False]
    k5 = f(t + dt, [u + 439*dt*k1[0]/216 - 8*dt*k2[0] + 3680*dt*k3[0]/513 - 845*dt*k4[0]/4104, v + 439*dt*k1[1]/216 - 8*dt*k2[1] + 3680*dt*k3[1]/513 - 845*dt*k4[1]/4104])
    
    if v - 8*dt*k1[1]/27 + 2*dt*k2[1] -3544*dt*k3[1]/2565 + 1859*dt*k4[1]/4104 -11*dt*k5[1]/40<0:
        return [0,0,0,0,False]
    k6 = f(t + dt/2, [u - 8*dt*k1[0]/27 + 2*dt*k2[0] -3544*dt*k3[0]/2565 + 1859*dt*k4[0]/4104 -11*dt*k5[0]/40, v - 8*dt*k1[1]/27 + 2*dt*k2[1] -3544*dt*k3[1]/2565 + 1859*dt*k4[1]/4104 -11*dt*k5[1]/40])

    du = 25*dt*k1[0]/216 + 1408*dt*k3[0]/2565 + 2197*dt*k4[0]/4101 - dt*k5[0]/5
    dv = 25*dt*k1[1]/216 + 1408*dt*k3[1]/2565 + 2197*dt*k4[1]/4101 - dt*k5[1]/5
    du1 = 16*dt*k1[0]/135 + 6656*dt*k3[0]/12825 + 28561*dt*k4[0]/56430 - 9*dt*k5[0]/50 + 2*dt*k6[0]/55
    dv1 = 16*dt*k1[1]/135 + 6656*dt*k3[1]/12825 + 28561*dt*k4[1]/56430 - 9*dt*k5[1]/50 + 2*dt*k6[1]/55
        
    return [du, dv,du1,dv1,True]

    
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
        
    return [y[0], y[1], dt, y[4]]


def ODEsolver(eq_type, central_pressure):
      
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
        
        if dy[3]==False:
            break
        mass = np.append(mass, mass[i] + dy[0])
        pressure = np.append(pressure, pressure[i] + dy[1])
        radius = np.append(radius, radius[i]+step)          
        i = i+1

        #print("step:",step,"radius:",radius[i],"mass:",mass[i],"pressure:",pressure[i])

    return radius,mass,pressure


def CubicSpline(x, y):
    
    #Split arrays into intervals
    n = len(x)
    h = []
    for i in range(n-1):
        h.append(x[i+1] - x[i])
        
    #To solve the system of equations of the spline we use the second derivative, that is linear in x for a cubic.
    #Integrating 2 times and using continuity we remove the 2n integration constants. Then, by derivating and matching
    #first derivatives we explicitate n-1 M. The remaining ones come from the boundary conditions that are set to 0.
    #The system to solve is made by a triangular matrix and we introduce c_p and d_p to solve it.

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
    
    #RHS of the system with tridiagonal
    D = [0] 
    for i in range(1, n-1):
        D.append(6*((y[i+1] - y[i])/h[i] - (y[i] - y[i-1])/h[i-1])/(h[i] + h[i-1]))
    D.append(0)
    
    #Solutions of the system with Thomas method
    c_p = C + [0] 
    d_p = [0]*n
    M = [0]*n 
    c_p[0] = C[0]/B[0]
    d_p[0] = D[0]/B[0]
    for i in range(1, n):
        c_p[i] = (c_p[i]/(B[i] - c_p[i-1]*A[i-1]))
        d_p[i] = (D[i] - d_p[i-1]*A[i-1])/(B[i] - c_p[i-1]*A[i-1])

    #M is the second derivative of the spline in all points, solution of the system above
    M[-1] = d_p[-1]
    for i in range(n-2, -1, -1):
        M[i] = d_p[i] - c_p[i]*M[i + 1]
    
    #Formula for cubic spline coefficients
    coefficients = []
    for i in range(n-1):
        coefficients.append([((M[i+1] - M[i])*h[i]**2)/6, (M[i]*h[i]**2)/2, (y[i+1] - y[i] - ((M[i+1] + 2*M[i])*h[i]**2)/6), y[i]])

    #bisect(list, num):This function returns the position in the sorted list, 
    #where the number passed in argument can be placed so as to maintain the resultant list in sorted order. 
    #If the element is already present in the list, the right most position where element has to be inserted is returned.
    def spline(val): 
        idx = min(bisect.bisect(x, val) - 1, n-1)
        z = (val - x[idx])/h[idx]
        Coef = coefficients[idx]
        S_z = Coef[0]*z**3 + Coef[1]*z**2 + Coef[2]*z + Coef[3]
        return S_z 
    return spline 


eos_library = {
#    'PAL6'  :[ 34.380,  2.227,  2.189,  2.159, 'npem' ],
    'SLy'   :[ 34.384,  3.005,  2.988,  2.851, 'npem' ],
#    'APR1'  :[ 33.943,  2.442,  3.256,  2.908, 'npem' ],
#    'APR2'  :[ 34.126,  2.643,  3.014,  2.945, 'npem' ],
    'APR3'  :[ 34.392,  3.166,  3.573,  3.281, 'npem' ],
    'APR4'  :[ 34.269,  2.830,  3.445,  3.348, 'npem' ],
#    'FPS'   :[ 34.283,  2.985,  2.863,  2.600, 'npem' ],
    'WFF1'  :[ 34.031,  2.519,  3.791,  3.660, 'npem' ],
    'WFF2'  :[ 34.233,  2.888,  3.475,  3.517, 'npem' ],
#    'WFF3'  :[ 34.283,  3.329,  2.952,  2.589, 'npem' ],
#    'BBB2'  :[ 34.331,  3.418,  2.835,  2.832, 'npem' ],
#    'BPAL12':[ 34.358,  2.209,  2.201,  2.176, 'npem' ],
    'ENG'   :[ 34.437,  3.514,  3.130,  3.168, 'npem' ],
    'MPA1'  :[ 34.495,  3.446,  3.572,  2.887, 'npem' ],
    'MS1'   :[ 34.858,  3.224,  3.033,  1.325, 'npem' ],
#    'MS2'   :[ 34.605,  2.447,  2.184,  1.855, 'npem' ],
    'MS1b'  :[ 34.855,  3.456,  3.011,  1.425, 'npem' ],
#    'PS'    :[ 34.671,  2.216,  1.640,  2.365, 'meson' ],
#    'GS1a'  :[ 34.504,  2.350,  1.267,  2.421, 'meson' ],
#    'GS2a'  :[ 34.642,  2.519,  1.571,  2.314, 'meson' ],
#    'BGN1H1':[ 34.623,  3.258,  1.472,  2.464, 'hyperon' ],
#    'GNH3'  :[ 34.648,  2.664,  2.194,  2.304, 'hyperon' ],
#    'H1'    :[ 34.564,  2.595,  1.845,  1.897, 'hyperon' ],
#    'H2'    :[ 34.617,  2.775,  1.855,  1.858, 'hyperon' ],
#    'H3'    :[ 34.646,  2.787,  1.951,  1.901, 'hyperon' ],
    'H4'    :[ 34.669,  2.909,  2.246,  2.144, 'hyperon' ],
#    'H5'    :[ 34.609,  2.793,  1.974,  1.915, 'hyperon' ],
#    'H6a'   :[ 34.593,  2.637,  2.121,  2.064, 'hyperon' ],
#    'H7'    :[ 34.559,  2.621,  2.048,  2.006, 'hyperon' ],
#    'PCL2'  :[ 34.507,  2.554,  1.880,  1.977, 'hyperon' ],
#    'ALF1'  :[ 34.055,  2.013,  3.389,  2.033, 'quark' ],
    'ALF2'  :[ 34.616,  4.070,  2.411,  1.890, 'quark' ],
#    'ALF3'  :[ 34.283,  2.883,  2.653,  1.952, 'quark' ],
#    'ALF4'  :[ 34.314,  3.009,  3.438,  1.803, 'quark' ]
}