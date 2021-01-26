import numpy as np
import scipy
import pyneb as pn
from ion_structure import compute_Ni,compute_NF,compute_NHIyi,compute_NHIy0
from cooling_rates import compute_lambda_CII_h,compute_lambda_CII_e

O3 = pn.Atom('O', 3) #load the OIII ion from pyneb

def oxygen_abundance(Z):
    """
    Set the oxygen abundance. 
    We assume Asplund et al 2009 abundance at Zsun and that Ao scales linearly with Z. Z in solar units
    """
    Ao = 4.90e-4
    return Ao*Z

def carbon_abundance(Z):
    """
    Set the carbon abundance. 
    We assume Asplund et al 2009 abundance at Zsun and that Ac scales linearly with Z. Z in solar units
    """
    Ac = 2.7e-4 
    return Ac*Z
   
def Sigmag_of_Sigmasfr(Sigma_sfr, k, ks_fit=1.4):
    """
    Gas surface density from the SFR surface density assuming the Kennicutt-Schmidt relation. 
    Sigma_sfr in Msun/yr/kpc^2, k is the burstiness parameter
    """
    out=(((k*(10.0**-12))**-1)*Sigma_sfr)**(1./ks_fit)
    return out

def U_Sigmag_Sigmasfr(Sigma_sfr, Sigma_g):
    """
    Ionizaton parameter as a function of star formation and gas surface densities. 
    reference:
      Eq. 38 Ferrara et al. 2019
    inputs:
      Sigma_sfr in Msun/yr/kpc^2
      Sigma_g   in Msun/kpc^2
    """
    out= (1.7e+14)*(Sigma_sfr/(Sigma_g*Sigma_g))
    return out

def compute_U_and_N(Z,k,Sigma_sfr,ks_fit=1.4 ):

    """
    wrapper to compute ionization parameter and column density
    """

    # the gas surface density is computed from the modified SK
    surface_density_gas = Sigmag_of_Sigmasfr(Sigma_sfr=Sigma_sfr, k=k,ks_fit=ks_fit)
    # convert the gas surface density to a column density
    column_density      = (surface_density_gas*10**22.0)/7.5e+7
    # compute the ionization parameter
    ionization_parameter= U_Sigmag_Sigmasfr(Sigma_sfr=Sigma_sfr, Sigma_g = surface_density_gas)

    return ionization_parameter,column_density 


def compute_flux_cii_density_bound(n, Z, U, column,  TPDR=100.0, THII=1.e+4):
    """
    Emerging [CII] flux for the density bounded case dens bounded case
    Eq. 34 in Ferrara et al. 2019
    """
    from atomic_data import g2_cii,g1_cii,E12_158um,A21_158um,ev2erg,n_crit_CII

    # column density in the density bound case
    N_HIy0 = compute_NHIy0(U, Z)

    column_neutral = column - N_i
    column_ionized = N_HIy0

    # compute the cooling contributions with and without LTE approximation
    if(n<=n_crit_CII):
        cooling_ionized    = n*compute_lambda_CII_e(T=THII)
    else:
        cooling_ionized    = (g2_cii/g1_cii)*pop_LTE(T_in=THII,E_in=E12_158um)*A21_158um*(ev2erg*E12_158um)

    out = carbon_abundance(Z)*(column_ionized*cooling_ionized)

    return out

def pop_LTE(T_in,E_in):

    """
    return the local thermodynamical equilibrium population
    convenient conversion for
      E_in is the level energy difference in eV
      T the temperature of the gas in K
    """

    from atomic_data import ev2erg

    out = np.exp((-ev2erg*E_in)/(1.38065e-16*T_in))

    return out

def compute_flux_cii_ionization_bound_N0(n, Z, U, column, TPDR=100.0, THII=1.e+4):

    """
    Emerging [CII] flux in the ionization bounded case, and for NF> N0
    Equations 30, 31 in Ferrara et al. 2019
    """

    from atomic_data import g2_cii,g1_cii,E12_158um,A21_158um,ev2erg,n_crit_CII

    # ionized column density
    N_i= compute_Ni(U,Z)
    # neutral column density in the ionized layer
    N_HIyi = compute_NHIyi(U, Z)

    column_neutral = column - N_i
    column_ionized = N_HIyi

    # compute the cooling contributions with and without LTE approximation
    if(n<=n_crit_CII):
        cooling_neutral    = n*compute_lambda_CII_h(T=TPDR)
        cooling_ionized    = n*compute_lambda_CII_e(T=THII)
    else:
        cooling_neutral    = (g2_cii/g1_cii)*pop_LTE(T_in=TPDR,E_in=E12_158um)*A21_158um*(ev2erg*E12_158um)
        cooling_ionized    = (g2_cii/g1_cii)*pop_LTE(T_in=THII,E_in=E12_158um)*A21_158um*(ev2erg*E12_158um)

    out = carbon_abundance(Z)*( column_neutral*cooling_neutral + column_ionized*cooling_ionized)

    return out

def compute_flux_cii_ionization_bound_NF(n, Z, U, column, TPDR=100.0, THII=1.e+4):
    """
    Emerging [CII] flux in the ionization bounded case, and for NF< N0
    Equations 30, 32 in Ferrara et al. 2019
    """

    from atomic_data import g2_cii,g1_cii,E12_158um,A21_158um,ev2erg,n_crit_CII

    # ionized column density
    N_i    = compute_Ni(U,Z)
    # NL column density
    N_F    = compute_NF(U, Z)
    # neutral column density in the ionized layer
    N_HIyi = compute_NHIyi(U, Z)

    column_neutral = N_F - N_i
    column_ionized = N_HIyi

    # compute the cooling contributions with and without LTE approximation
    if(n<=n_crit_CII):
        cooling_neutral    = n*compute_lambda_CII_h(T=TPDR)
        cooling_ionized    = n*compute_lambda_CII_e(T=THII)
    else:
        cooling_neutral    = (g2_cii/g1_cii)*pop_LTE(T_in=TPDR,E_in=E12_158um)*A21_158um*(ev2erg*E12_158um)
        cooling_ionized    = (g2_cii/g1_cii)*pop_LTE(T_in=THII,E_in=E12_158um)*A21_158um*(ev2erg*E12_158um)

    out = carbon_abundance(Z)*( column_neutral*cooling_neutral + column_ionized*cooling_ionized)
        
    return out

def Sigma_CII158(logn, Z, k, Sigma_sfr):

    """
    compute the CII surface brightenss
    Eq. 35 in Ferrara et al. 2019
    """

    UU , column_density = compute_U_and_N(Z=Z,k=k,Sigma_sfr=Sigma_sfr)
    N_i                 = compute_Ni(U=UU,Z=Z)
    N_F                 = compute_NF(U=UU,Z=Z)
    n                   = 10**logn

    if(column_density<N_i):
         flux = compute_flux_cii_density_bound(n = n, Z= Z, U = UU, column=column_density   , TPDR=100.0, THII=1e+4)
    elif(column_density<N_F):
         flux = compute_flux_cii_ionization_bound_N0(n = n, Z= Z, U = UU, column=column_density, TPDR=100.0, THII=1e+4)
    else:
         flux = compute_flux_cii_ionization_bound_NF(n = n, Z= Z, U = UU, column=column_density, TPDR=100.0, THII=1e+4)

    # flux to surface brightness conversion
    out = flux*2.474e+9

    return out

def foiii(n, Z, U, THII,line="88um"):

    """
    flux for the [OIII] line emission (88 and 52 micron), 
    Emissivity computed with Pyneb.
    details can be found in Vallini et al. 2021.
    """

    emOIII = None
    # get emissivity and relative population for the selected line:  
    if line == "52um":
      emOIII        = O3.getEmissivity(THII, n, wave='51.8m') # erg s^-1 cm^3
    if line == "88um":
      emOIII        = O3.getEmissivity(THII, n, wave='88.3m') # erg s^-1 cm^3

    # ionized column density
    N_i           = compute_Ni(U=U,Z=Z)
    #correction for the presence of OII in the ionized region
    fo3           = np.array([0.10994503, 0.73298314, 0.96966708])
    Uo3           = np.array([-3.5, -2.5, -1.5])
    Xoiii         = np.interp(np.log10(U), Uo3, fo3)

    Nh_oiii       = oxygen_abundance(Z)*Xoiii*N_i
    out           = emOIII * n * Nh_oiii

    return out

def Sigma_OIII88(logn, Z, k, Sigma_sfr):
    UU , __ =  compute_U_and_N(Z=Z,k=k,Sigma_sfr=Sigma_sfr)
    ff      = foiii(n=10.0**logn, Z=Z, U=UU, THII=1e+4,line="88um")
    out     = ff*2.474e+9 # 
    return out

def Sigma_OIII52(logn, Z, k, Sigma_sfr):
    UU , __ =  compute_U_and_N(Z=Z,k=k,Sigma_sfr=Sigma_sfr)
    ff      = foiii(n=10.0**logn, Z=Z, U=UU, THII=1e+4,line="52um")
    out     = ff*2.474e+9
    return out

def Delta(logn, Z, k, Sigma_sfr):

    from empirical import delooze_fit_resolved
    
    out = np.log10(Sigma_CII158(logn=logn, Z=Z, k=k, Sigma_sfr=Sigma_sfr))-np.log10(delooze_fit_resolved(Sigma_sfr=Sigma_sfr))
    return out




