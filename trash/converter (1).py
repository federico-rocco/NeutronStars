import numpy as np 

c = 3*10**(10)

G = 6.67*10**(-7)

M_sun = 2*10**(33)

def from_cgs_to_geo(m, e, ed, md, p):
    m = m*G/(c**2)
    md = md*G/(c**2)
    e = e*G/(c**4)
    ed = ed*G/(c**4)
    p = p*G/(c**4)
    return m, md, e, ed, p
    
def from_geo_to_cgs(m, e, ed, md, p):
    m = m*(c**2)/G
    md = md*(c**2)/G
    e = e*(c**4)/G
    ed = ed*(c**4)/G 
    p = p*(c**4)/G
    return m, md, e, ed, p

def mass_conv(m):
    m = m/M_sun
    return m

def radius_conv(r):
    r = r*10**(-6)
    











