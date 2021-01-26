import numpy as np

"""
    this module contains the analytical eq.s for the cooling rate of various ions
    outpus:
      cooling rates in erg/s/cm-3
    inputs:
      T, the temperature in K
"""

def lambdaCIIh(T):
    
    """
    Maxwellian-averaged collision rates with neutrals
    using expression from Goldsmith et al. 2012
    see Appendix A, Ferrara et al. 2019
    """

    factor = (1.84e-4*(T)**0.64)/2.0
    out    = (8.6293e-6/np.sqrt(T))*factor*(1.602e-12*0.0079)*np.exp((-1.602e-12*0.0079)/(1.38065e-16*T))

    return out

def lambdaCIIe(T):
    """
    Maxwellian-averaged collision rates with e-
    using expression from Goldsmith et al. 2012.
    see Appendix A, Ferrara et al. 2019
    """
    factor= (0.67*(T)**0.13)/2.0
    out   = (8.6293e-6/np.sqrt(T))*factor*(1.602e-12*0.0079)*np.exp((-1.602e-12*0.0079)/(1.38065e-16*T))

    return out



