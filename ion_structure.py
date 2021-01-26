
import numpy as np

"""
  this module contains the equations to compute the ionization structure
  inputs:
    U is the ionzation parameter
    Z is the metallicity in Zsun
  outputs:
    column densities, measured in cm^-2
"""

def Ni(U,Z):
    """
    Column density of the ionized gas
    Eq. 14 in Ferrara et al. 2019. 
    """
    ND=Nd(Z)
    tau_sd = (1e+23*U)/ND
    NN=ND*np.log((1+tau_sd)/(1+(tau_sd/np.exp(1.0))))
    return NN

def Nd(Z):
    """
    Column density corresponding to A_V = 1 
    Eq. 9c in Ferrara et al. 2019
    """
    Ndsolar=1.7e21 #cm^-2
    return Ndsolar/Z

def chi_of_U(Z):
    """
    Eq. 22 in Ferrara et al. 2019
    """
    
    chi = 8.7e+4*Z

    return chi

def chi_prime(U, Z):
    """
    Eq. 22b in Ferrara et al. 2019
    """
    ww       = w_of_D(Z)
    chi      = chi_of_U(U)
    chiprime = ww*chi
    return chiprime

def w_of_D(Z):
    """
    Factor related to the abs. of LW photons
    see Eq. 24 in Ferrara et al. 2019 (see also Sternberg et al. 2014)
    """
    w = 1.0/(1.0+ 0.9*(Z)**0.5)
    return w

def NF(U,Z):
    """
    Column density at which the Lyman-Werner flux vanishes
    Eq. 28 Ferrara et al. 2019
    """
    ND = Nd(Z)
    out = ND*np.log(1+chi_prime(U,Z))
    return out

def NHIyi(U, Z):
    """
    HI column in the ionized layer
    Eq. 13 and 14 in Ferrara et al. 2019
    """
    Ns     = 1e+23*U
    tau_s  = 2.7e+5 * U
    out    = (Ns/tau_s) * (1.0 - Ni(U,Z)/Nd(Z))
    return out

def NHIy0(U, Z, N0):
    """
    Ionized column density in the case of the density bounded regime
    Eq. 33 in Ferrara et al. 2019
    """
    Ns      = 1e+23*U
    y0      = N0/Ns
    tau_sd  = 59.0 * U * Z
    tau_s   = 2.7e+5 * U
    out     = (Ns/tau_s) * np.log(tau_sd/np.abs((np.exp(tau_sd*y0)-tau_sd-1.0)))
    return out


