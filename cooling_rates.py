import numpy as np


# Maxwellian-averaged collision rates with neutrals (Appendix A, Ferrara et al. 2019) 
#using expression from Goldsmith et al. 2012. Temperature in Kelvin
def lambdaCIIh(T):
    factor=(1.84e-4*(T)**0.64)/2.0
    out = (8.6293e-6/np.sqrt(T))*factor*(1.602e-12*0.0079)*np.exp((-1.602e-12*0.0079)/(1.38065e-16*T))
    return out

# Maxwellian-averaged collision rates with e- (Appendix A, Ferrara et al. 2019) 
#using expression from Goldsmith et al. 2012. Temperature in Kelvin
def lambdaCIIe(T):
    factor=(0.67*(T)**0.13)/2.0
    out = (8.6293e-6/np.sqrt(T))*factor*(1.602e-12*0.0079)*np.exp((-1.602e-12*0.0079)/(1.38065e-16*T))
    return out



