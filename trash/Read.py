# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 11:04:59 2021

@author: lucab
"""

#prova classe Read 

"""
COSA FARE A PAROLE: 

il modello di neutron star è diviso fondamentalmente in due parti:
    
    OUTER REGION: basata su SLY, uguale per tutti e sempre con gli stessi parametri
    
    CORE REGION: sono i modelli che devo creare a manetta, ognuno di essi ha:
        
        tre esponenti, un valore di matching, un errore da riportare 
        
        tre coefficienti da calcolare
        
        due intersezioni da capire come determinare
        
        una densità centrale che deve variare in un certo intervallo 

COME FARE IN PRATICA:
    
creo una classe Read

le caratteristiche degli oggetti sono quelli che devo dargli io fuori dalla classe

metto dentro la classe lo SLY per la inner region (riguardare l'esempio numero ore)
                                                   
metto anche un metodo che agisce sull'oggetto per i tre coefficienti 

metto un altro metodo per l'eventuale intersezione da determinare

un parametro dovrà essere libero, la densità centrale da buttare in un ciclo for...?

ESEMPIO: 

Generic polytropic essere tipo p = K*rho**Gamma

PAL6:

log(p1) = 34.380 è la rho_in/out

Gamma_in_1 = 2.227

Gamma_in_2 = 2.189

Gamma_in_3 = 2.159

K_in_1 = p(rho_in/out)/rho_in/out**(Gamma_out_3)

K_in_2 = p(rho_in_1)/rho_in_1**(Gamma_in_1)

K_in_3 = p(rho_in_2)/rho_in_2**(Gamma_in_2)

rho_in_1 e rho_in_2 sono date nell'articolo
central density a parte va preso da compact stars
"""
class stella_di_read: 
    
    #outer region, SLY model
    K_sly_0 = 6.801*10**(-9)
    K_sly_1 = 1.061*10**(-6)
    K_sly_2 = 5.326*10**(1)
    K_sly_3 = 3.998*10**(-8)
    gamma_sly_0 = 1.584
    gamma_sly_1 = 1.287
    gamma_sly_2 = 0.622
    gamma_sly_3 = 1.356
    rho_sly_0 = 0
    rho_sly_1 = 2.440*10**7
    rho_sly_2 = 3.783*10**11
    rho_sly_3 = 2.627*10**12 
    
    #interior region, particular model
    def __init__(self, gamma_in_1, gamma_in_2, gamma_in_3, p1, central_density):
        self.gamma_in_1 = gamma_in_1 
        self.gamma_in_2 = gamma_in_2
        self.gamma_in_3 = gamma_in_3
        self.p1 = p1
        self.central_density = central_density
    
    def k_calculation(prev_dens, prev_pres, resp_gamma):
        return prev_pres/((prev_dens)**(resp_gamma))
    
    def pres_poly(rho, K, Gamma):
        return K*rho**(Gamma)
    rho_in/out = 1 ##per ora gli do un valore, ma va tirata fuori da p1
    #non penso siano queste ma per ora mettiamole, devono essere fixate per tutti....
    rho_in_2 = 10**(14)
    rho_in_3 = 10**(15)
 ##come chiamo la variabile?   
    def equation_of_state(rho):
        if rho_sly_0 < rho < rho_sly_1 :
            return K_sly_0*rho**(gamma_sly_0)
        elif rho_sly_1 < rho < rho_sly_2 :
            return K_sly_1*rho**(gamma_sly_1)
        elif rho_sly_2 < rho < rho_sly_3:
            return K_sly_2*rho**(gamma_sly_2)
        elif rho_sly_3 < rho < rho_in/out:
            return K_sly_3*rho**(gamma_sly_3)
        elif rho_in/out < rho < rho_in_2:
            return k_calculation()*rho**(gamma_in_1)
        elif rho_in_2 < rho < rho_in_3:
            return k_calculation(rho_in_2, pres_poly(rho_in_2, k_calculation(AIUTOOO), gamma_in_2))*rho**(gamma_in_2)
        else: 
            return k_calculation(rho_in_3, pres_poly(rho_in_3, k_calculation(), gamma_in_3))*rho**(gamma_in_3)
#attento, potrebbe essere sbagliata la rho nel k calculation, non sono sicuro se è 2 o 1 per la prima e 3 o 2 per la seconda        
        
#AAAAAAAAAAA TROPPE FUNZIONI NELLE FUNZIONI NELLE FUNZIONI NELLE FUNZIONI NON SO PROGRAMMARE AIUTO FEDE
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        