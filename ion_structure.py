
import numpy as np

# Column density of the ionized gas, Eq. 14 in Ferrara et al. 2019. 
#U is the ion. parameter, Z is the metallicity in Zsun.
def Ni(U,Z):
    ND=Nd(Z)
    tau_sd = (1e+23*U)/ND
    NN=ND*np.log((1+tau_sd)/(1+(tau_sd/np.exp(1.0))))
    return NN

# Column density corresponding to A_V = 1, as a function of Z. 
#Eq. 9c in Ferrara et al. 2019
def Nd(Z):
    Ndsolar=1.7e21 #cm^-2
    return Ndsolar/Z

# Eq. 22 in Ferrara et al. 2019
def chi_of_U(Z):
    chi=8.7e+4*Z
    return chi

# Eq. 22b in Ferrara et al. 2019
def chi_prime(U, Z):
    ww=w_of_D(Z)
    chi=chi_of_U(U)
    chiprime=ww*chi
    return chiprime

# Factor related to the abs. of LW photons, see Eq. 24 in Ferrara et al. 2019 (see also Sternberg et al. 2014)
def w_of_D(Z):
    w=1.0/(1.0+ 0.9*(Z)**0.5)
    return w

# Column density at which the Lyman-Werner flux vanishes, Eq. 28 Ferrara et al. 2019
def NF(U,Z):
    ND=Nd(Z)
    NN=ND*np.log(1+chi_prime(U,Z))
    return NN

# HI column in the ionized layer. Eq. 13 and 14 in Ferrara et al. 2019
def NHIyi(U, Z):
    Ns=1e+23*U
    tau_s  = 2.7e+5 * U
    out = (Ns/tau_s) * (1.0 - Ni(U,Z)/Nd(Z))
    return out

# Eq. 33 in Ferrara et al. 2019. Ionized column density in the case of the density bounded regime
def NHIy0(U, Z, N0):
    Ns=1e+23*U
    y0=N0/Ns
    tau_sd = 59.0 * U * Z
    tau_s  = 2.7e+5 * U
    out = (Ns/tau_s) * np.log(tau_sd/np.abs((np.exp(tau_sd*y0)-tau_sd-1.0)))
    return out


