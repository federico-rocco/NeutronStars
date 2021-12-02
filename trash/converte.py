import numpy as np 

c = 3*10**(10)

G = 6.67*10**(-7)



def from_cgs_to_geo(quantity, value):
    if quantity == "mass":
        return value*G/(c**2)
    elif quantity == "massden":
        return value*G/(c**2)
    elif quantity == "energy":
        return value*G/(c**4)
    elif quantity == "eden":
        return value*G/(c**4)
    elif quantity == "pressure":
        return value*G/(c**4)
    else:
        raise ValueError("Need to choose \"mass\" or \"massden\" or \"energy\" or \"eden\" or \"pressure\"")
        
def from_geo_to_cgs(quantity, value):
    if quantity == "mass":
        return value*(c**2)/G
    elif quantity == "massden":
        return value*(c**2)/G
    elif quantity == "energy":
        return value*(c**4)/G
    elif quantity == "eden":
        return value*(c**4)/G 
    elif quantity == "pressure":
        return value*(c**4)/G
    else:
        raise ValueError("Need to choose \"mass\" or \"massden\" or \"energy\" or \"eden\" or \"pressure\"")

def mass_conv(m):
    M_sun = 2*10**(33)
    m = m/M_sun
    return m

def radius_conv(r):
    r = r*10**(-6)
    











