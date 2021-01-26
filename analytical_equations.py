import numpy as np
import scipy
import pyneb as pn

O3 = pn.Atom('O', 3) #load the OIII ion from pyneb

# Set the oxygen abundance. 
#We assume Asplund et al 2009 abundance at Zsun and that Ao scales linearly with Z. Z in solar units
def oxygen_abundance(Z):
    Ao=4.90e-4
    return Ao*Z

# Set the carbon abundance. 
#We assume Asplund et al 2009 abundance at Zsun and that Ac scales linearly with Z. Z in solar units
def carbon_abundance(Z):
    Ac=2.7e-4 
    return Ac*Z
     
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

# Gas surface density from the SFR surface density assuming the Kennicutt-Schmidt relation. 
#Sigma_sfr in Msun/yr/kpc^2, k is the burstiness parameter
def Sigmag_of_Sigmasfr(Sigma_sfr, k, n=1.4):
    out=(((k*(10.0**-12))**-1)*Sigma_sfr)**(1./n)
    return out

# Ionizaton parameter from SFR and gas surface densities. 
#Eq. 38 Ferrara et al. 2019. Sigma_sfr in Msun/yr/kpc^2, Sigma_g in Msun/kpc^2
def U_Sigmag_Sigmasfr(Sigma_sfr, Sigma_g):
    out= (1.7e+14)*(Sigma_sfr/(Sigma_g*Sigma_g))
    return out

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

# Emerging [CII] flux for the density bounded case dens bounded case.  Eq. 34 in Ferrara et al. 2019
def fcii_DB(n, Z, U, column, TPDR, THII):
    g2_cii=4. # 
    g1_cii=2.
    E12_158um = 0.0079 #energy between the eV
    A21_158um = 2.4e-6 #s^-1
    if(n<=3300):
        # rates:  
        lambdaCII    = lambdaCIIh(TPDR)
        lambdaCII4   = lambdaCIIe(THII)

        fcii_neutral = 0.0
        fcii_ionized_DB=n*carbon_abundance(Z)*lambdaCII4*NHIy0(U, Z, column)
    
        out= fcii_neutral + fcii_ionized_DB
    else:
        
        fcii_neutral = 0.0
        LTE_pop_levels_HII = (g2_cii/g1_cii)*np.exp((-1.602e-12*E12_158um)/(1.38065e-16*THII))
        fcii_ionized_DB = LTE_pop_levels_HII * carbon_abundance(Z) * A21_158um * (1.602e-12*E12_158um)* NHIy0(U, Z, column)
        out = fcii_neutral + fcii_ionized_DB
        
    return out

# Emerging [CII] flux in the ionization bounded case, and for NF> N0. Equations 30, 31 in Ferrara et al. 2019
def fcii_IB_N0(n, Z, U, column, TPDR, THII):
    g2_cii=4.
    g1_cii=2.
    E12_158um = 0.0079 #eV
    A21_158um = 2.4e-6 #s^-1
    # ionized column density
    N_i= Ni(U,Z)
    
    if(n<=3300):
        # rates:  
        lambdaCII    = lambdaCIIh(TPDR)
        lambdaCII4   = lambdaCIIe(THII)

        # cii from neutral
        fcii_neutral = n*carbon_abundance(Z)*lambdaCII*(column - N_i)
    
        # cii from ionized layer
        fcii_ionized_IB=n*carbon_abundance(Z)*lambdaCII4*NHIyi(U, Z)
    
        out= fcii_neutral + fcii_ionized_IB
    else:
        LTE_pop_levels_PDR = (g2_cii/g1_cii)*np.exp((-1.602e-12*E12_158um)/(1.38065e-16*TPDR))
        LTE_pop_levels_HII = (g2_cii/g1_cii)*np.exp((-1.602e-12*E12_158um)/(1.38065e-16*THII))
        fcii_neutral = LTE_pop_levels_PDR * carbon_abundance(Z) * A21_158um * (1.602e-12*E12_158um)* (column - N_i)
        fcii_ionized_IB = LTE_pop_levels_HII * carbon_abundance(Z) * A21_158um * (1.602e-12*E12_158um)* NHIyi(U, Z)
        out= fcii_neutral + fcii_ionized_IB
        
    return out

# Emerging [CII] flux in the ionization bounded case, and for NF< N0. Equations 30, 32 in Ferrara et al. 2019
def fcii_IB_NF(n, Z, U, column, TPDR, THII):
    g2_cii=4.
    g1_cii=2.
    E12_158um = 0.0079 #eV
    A21_158um = 2.4e-6 #s^-1
    # ionized column density
    N_i= Ni(U,Z)
    # NL column density
    N_F=NF(U, Z)
    
    if(n<=3300):
        # rates:  
        lambdaCII    = lambdaCIIh(TPDR)
        lambdaCII4   = lambdaCIIe(THII)

        # cii from neutral
        fcii_neutral = n*carbon_abundance(Z)*lambdaCII*(N_F - N_i)
    
        # cii from ionized layer
        fcii_ionized_IB=n*carbon_abundance(Z)*lambdaCII4*NHIyi(U, Z)
    
        out= fcii_neutral + fcii_ionized_IB
        
    else:
        
        LTE_pop_levels_PDR = (g2_cii/g1_cii)*np.exp((-1.602e-12*E12_158um)/(1.38065e-16*TPDR))
        LTE_pop_levels_HII = (g2_cii/g1_cii)*np.exp((-1.602e-12*E12_158um)/(1.38065e-16*THII))
        fcii_neutral = LTE_pop_levels_PDR * carbon_abundance(Z) * A21_158um * (1.602e-12*E12_158um)* (N_F - N_i)
        fcii_ionized_IB = LTE_pop_levels_HII * carbon_abundance(Z) * A21_158um * (1.602e-12*E12_158um)* NHIyi(U, Z)
        out= fcii_neutral + fcii_ionized_IB
        
    return out


# The following three equations are instrumental for the Eq. 35 in Ferrara et al. 2019.
def sigma_cii_DB(logn, Z, k, Sigma_sfr):
    n= 10**logn
    SS_g=Sigmag_of_Sigmasfr(Sigma_sfr, k)
    UU=U_Sigmag_Sigmasfr(Sigma_sfr, SS_g)
    column_density=(SS_g*10**22.0)/7.5e+7
    ff=fcii_DB(n, Z, UU, column_density, 100., 1e+4)
    SS_CII=ff*2.474e+9
    return SS_CII

def sigma_cii_IB_N0(logn, Z, k, Sigma_sfr):
    n= 10**logn
    SS_g=Sigmag_of_Sigmasfr(Sigma_sfr, k)
    UU=U_Sigmag_Sigmasfr(Sigma_sfr, SS_g)
    column_density=(SS_g*10**22.0)/7.5e+7
    ff=fcii_IB_N0(n, Z, UU, column_density, 100., 1e+4)
    SS_CII=ff*2.474e+9
    return SS_CII

def sigma_cii_IB_NF(logn, Z, k, Sigma_sfr):
    n= 10**logn
    SS_g=Sigmag_of_Sigmasfr(Sigma_sfr, k)
    UU=U_Sigmag_Sigmasfr(Sigma_sfr, SS_g)
    column_density=(SS_g*10**22.0)/7.5e+7
    ff=fcii_IB_NF(n, Z, UU, column_density, 100., 1e+4)
    SS_CII=ff*2.474e+9
    return SS_CII

# Eq. 35 in Ferrara et al. 2019
def Sigma_CII158(logn, Z, k, Sigma_sfr):
    
    SS_g=Sigmag_of_Sigmasfr(Sigma_sfr, k)
    UU=U_Sigmag_Sigmasfr(Sigma_sfr, SS_g)
    column_density=(SS_g*10**22.0)/7.5e+7
    N_i=Ni(UU,Z)
    N_F=NF(UU,Z)
    
    if(column_density<N_i):
        out = sigma_cii_DB(logn, Z, k, Sigma_sfr)
    elif(column_density<N_F):
            out=sigma_cii_IB_N0(logn, Z, k, Sigma_sfr)
    else:
            out=sigma_cii_IB_NF(logn, Z, k, Sigma_sfr)
    return out
    
# Part related to the [OIII] line emission (88 and 52 micron), 
#details can be found in Vallini et al. 2021. Emissivity computed with Pyneb.
def foiii88(n, Z, U, THII):
    # emissivity for 88micron:  
    emOIII  = O3.getEmissivity(THII, n, wave='88.3m') # erg s^-1 cm^3
    # ionized column density
    N_i=Ni(U,Z)
    
    #correction for the presence of OII in the ionized region
    fo3=np.array([0.10994503, 0.73298314, 0.96966708])
    Uo3=np.array([-3.5, -2.5, -1.5])
    Xoiii=np.interp(np.log10(U), Uo3, fo3)
    n1_ntot= O3.getPopulations(THII, n)[1]
    Nh_oiii = oxygen_abundance(Z)*Xoiii*N_i
    foiii_ionized = emOIII * n * Nh_oiii
    return foiii_ionized

def foiii52(n, Z, U, THII):
    #emissivity for 52micron:  
    emOIII  = O3.getEmissivity(THII, n, wave='51.8m') #erg s^-1 cm^3
    # ionized column density
    N_i=Ni(U,Z)
    #correction for the presence of OII in the ionized region
    fo3=np.array([0.10994503, 0.73298314, 0.96966708])
    Uo3=np.array([-3.5, -2.5, -1.5])
    Xoiii=np.interp(np.log10(U), Uo3, fo3)
    n2_ntot= O3.getPopulations(THII, n)[2]
    Nh_oiii = oxygen_abundance(Z)*Xoiii*N_i
    foiii_ionized = emOIII * n * Nh_oiii
    return foiii_ionized 

def Sigma_OIII88(logn, Z, k, Sigma_sfr):
    n=10**logn
    SS_g=Sigmag_of_Sigmasfr(Sigma_sfr, k)
    UU=U_Sigmag_Sigmasfr(Sigma_sfr, SS_g)
    ff=foiii88(n, Z, UU, 1e+4)
    out=ff*2.474e+9 # 
    return out

def Sigma_OIII52(logn, Z, k, Sigma_sfr):
    n=10**logn
    SS_g=Sigmag_of_Sigmasfr(Sigma_sfr, k)
    UU=U_Sigmag_Sigmasfr(Sigma_sfr, SS_g)
    ff=foiii52(n, Z, UU, 1e+4)
    out=ff*2.474e+9
    return out

def delooze_fit(logsfr):
    return (logsfr + 7.06)/1.0

def delooze_fit_resolved(Sigma_sfr):
    logSigma_cii=(np.log10(Sigma_sfr) +6.99)/0.93
    return 10**logSigma_cii

def delooze_delta(Sigma_sfr,Sigma_cii):
    return np.log10(Sigma_cii) - np.log10(delooze_fit_resolved(Sigma_sfr))

def Delta(logn, Z, k, Sigma_sfr):
    out = np.log10(Sigma_CII158(logn, Z, k, Sigma_sfr))-np.log10(delooze_fit_resolved(Sigma_sfr))
    return out




