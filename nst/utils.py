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



#dictionary for units provided by LIGO https://github.com/lscsoft/bilby/blob/master/bilby/gw/eos/eos.py
cgs_geom_dictionary = { "geom": { "lenght": {"cm": 100.,
                                             "m": 1,
                                             "km": 1e-3
                                             },
                                 
                                  "time": {"cgs": 1/C_SI,
                                           "si": 1/C_SI,
                                           "geom": 1.
                                           },
                                 
                                  "mass": { "g": 1e3*(C_SI**2.)/G_SI,
                                            "kg": (C_SI**2.)/G_SI,
                                            "m_sol": (C_SI**2.)/G_SI/MSUN_SI
                                           },
                                                                   
                                  "energy": {"cgs": 1e7*(C_SI**4.)/G_SI, 
                                             "si": (C_SI**4.)/G_SI,
                                             },
                                  
                                  "energy_density": {"cgs": 10*(C_SI**4.)/G_SI
                                                     },
                                  
                                  "density": {"cgs": 1e-3*(C_SI**2.)/G_SI
                                              },
                                  
                                  "pressure": {"cgs": 10*(C_SI**4.)/G_SI
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
                                          "geom": 1e-3*G_SI/(C_SI**2.) ,
                                          "m_sol": 1/MSUN_CGS
                                          },
                                 
                                  "energy": {"si": 1e-7, 
                                            "geom": 1e-7*G_SI/(C_SI**4.)
                                            },
                                                         
                                 "energy_density": {"si": 0.1,
                                                    "geom": 0.1*G_SI/(C_SI**4.) 
                                                    },
                                 
                                 "density": {"si": 1e3,
                                             "geom": 1e3*G_SI/(C_SI**2.)
                                             },
                                                         
                                 "pressure": {"si": 0.1,
                                              "geom": 0.1*G_SI/(C_SI**4.)
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
                                          "geom":  G_SI/(C_SI**2.),
                                          "m_sol": 1/MSUN_SI
                                          },
                                                                                                
                                 "energy": {"cgs": 1e7, 
                                            "geom":  G_SI/(C_SI**4.)
                                           },
                                
                                "energy_density": {"cgs": 10.,
                                                   "geom": G_SI/(C_SI**4.)
                                                   },
                                
                                 "density": {"cgs": 1e-3,
                                             "geom": (C_SI**2.)/G_SI,
                                             },
                                
                                "pressure": {"cgs": 10.,
                                             "geom": G_SI/(C_SI**4.)
                        
                                             },
                                }
                        }


    
def runge_kutta(f, t, u, v, dt):

    """
    
    Function that implements the 5th order Runge-Kutta algorithm

    Parameters
    ----------
    f : function
        system of two differential equations (in our case : TOV or Newton equations)
    t : float
        independent variable (in our case: radius)
    u : float
        previous value of the first dependent variable (in our case: mass)
    v : float
        previous value of the second dependent variable (in our case: pressure)
    dt : float
        step

    Returns [du, dv] : array of float
        array containing the increments of the dependent variables
    -------

    """    
    
    k1 = f(t, [u, v])
    k2 = f(t + dt/4, [u + dt*k1[0]/4, v + dt*k1[1]/4])
    k3 = f(t + 3*dt/8, [u + 3*dt*k1[0]/32 + 9*dt*k2[0]/32, v + 3*dt*k1[1]/32 + 9*dt*k2[1]/32])
    k4 = f(t + 12*dt/13, [u + 1932*dt*k1[0]/2197 - 7200*dt*k2[0]/2197 + 7296*dt*k3[0]/2197, v + 1932*dt*k1[1]/2197 - 7200*dt*k2[1]/2197 + 7296*dt*k3[1]/2197])
    k5 = f(t + dt, [u + 439*dt*k1[0]/216 - 8*dt*k2[0] + 3680*dt*k3[0]/513 - 845*dt*k4[0]/4104, v + 439*dt*k1[1]/216 - 8*dt*k2[1] + 3680*dt*k3[1]/513 - 845*dt*k4[1]/4104])

    du = 25*dt*k1[0]/216 + 1408*dt*k3[0]/2565 + 2197*dt*k4[0]/4101 - dt*k5[0]/5
    dv = 25*dt*k1[1]/216 + 1408*dt*k3[1]/2565 + 2197*dt*k4[1]/4101 - dt*k5[1]/5
                        
    return [du, dv]



def ode_solver(eq_type, initial_value, step=10):

    """
    
    Function that uses the Runge-Kutta algorithm with fixed step to resolve a system of differential equations
    A variable step approach didn't provide improved results, despite complicating the code
    Cauchy conditions:
        y_2(x=0) = initial_value
        y_2(x_tot) = 0

    Parameters
    ----------
    eq_type : function
        the system of differential equations to be solved
    initial_value : float
        Cauchy condition for y_2 at x = 0
    step : float, optional
        fixed step chosen for the computaion. The default is 10

    Returns x, y_1, y_2 : numpy arrays
        arrays containing the values of the independent and dependent variables at each step
    -------
    None.

    """  
    
    x = np.array([1e-10])
    y_1 = np.array([0])
    y_2 = np.array([initial_value])
    
    
    i=0
    
    while(y_2[i]>0):

        dy = runge_kutta(eq_type, x[i], y_1[i], y_2[i], step)
        y_1 = np.append(y_1, y_1[i] + dy[0])
        y_2 = np.append(y_2, y_2[i] + dy[1])
        x = np.append(x, x[i] + step)
      
        i += 1

    return x, y_1, y_2



def cubic_spline(x, y):

    """
    
    Function that, given two array of values of x and y computed in the same points, computes the interpolated solution 
    y(x) by implementing the cubic_spline algorithm

    Parameters
    ----------
    x : array of float
        array containing the values of the first function
    y : array of float
        array containing the values of the first function

    Returns spline : function
        y(x)
    -------

    """
    
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
    'SLy'   :[ 34.384,  3.005,  2.988,  2.851, 'npem' ],
    'APR3'  :[ 34.392,  3.166,  3.573,  3.281, 'npem' ],
    'APR4'  :[ 34.269,  2.830,  3.445,  3.348, 'npem' ],
    'WFF1'  :[ 34.031,  2.519,  3.791,  3.660, 'npem' ],
    'WFF2'  :[ 34.233,  2.888,  3.475,  3.517, 'npem' ],
    'ENG'   :[ 34.437,  3.514,  3.130,  3.168, 'npem' ],
    'MPA1'  :[ 34.495,  3.446,  3.572,  2.887, 'npem' ],
    'MS1'   :[ 34.858,  3.224,  3.033,  1.325, 'npem' ],
    'MS1b'  :[ 34.855,  3.456,  3.011,  1.425, 'npem' ],
    'H4'    :[ 34.669,  2.909,  2.246,  2.144, 'hyperon' ],
    'ALF2'  :[ 34.616,  4.070,  2.411,  1.890, 'quark' ]
}