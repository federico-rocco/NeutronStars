# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 11:19:33 2021

@author: lucab
"""
#https://en.wikiversity.org/wiki/Cubic_Spline_Interpolation

import bisect

def calculate_spline(x, y):
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