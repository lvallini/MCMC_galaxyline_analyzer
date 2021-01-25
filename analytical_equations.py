import numpy as np
import scipy
import pyneb as pn

O3 = pn.Atom('O', 3) #load the OIII ion from pyneb
C2 = pn.Atom('C', 2) #load the CII ion from pyneb

# oxygen abundance
def oxygen_abundance(Z):
    Ao=4.90e-4 #oxy abund at solar metallicity
    return Ao*Z

# carbon abundance
def carbon_abundance(Z):
    Ac=2.7e-4 #carbon abund at solar metallicity
    return Ac*Z
      
# collisional excitation CIII
def lambdaCIIIe(T):
    factor=1.265-0.0174e-4*T
    out = (8.6293e-6/np.sqrt(T))*factor*(1.602e-12*6.54)*np.exp((-1.602e-12*6.54)/(1.38065e-16*T))
    return out

# collisional excitation CII by H
def lambdaCIIh(T):
    factor=(1.84e-4*(T)**0.64)/2.0
    out = (8.6293e-6/np.sqrt(T))*factor*(1.602e-12*0.0079)*np.exp((-1.602e-12*0.0079)/(1.38065e-16*T))
    return out

# collisional excitation CII by e-
def lambdaCIIe(T):
    factor=(0.67*(T)**0.13)/2.0
    out = (8.6293e-6/np.sqrt(T))*factor*(1.602e-12*0.0079)*np.exp((-1.602e-12*0.0079)/(1.38065e-16*T))
    return out

# gas surface density from the SFR surface density 
def Sigmag_of_Sigmasfr(Sigma_sfr, k, n=1.4):
    out=(((k*(10.0**-12))**-1)*Sigma_sfr)**(1./n)
    return out

# gas surface density from the SFR surface density 
def Sigmag_of_Sigmasfr2(Sigma_sfr, k, n=1.4):
    A=k*(1.6e-4)
    Sigmagas_=(Sigma_sfr/A)**(1./n)
    out = Sigmagas_ * 1e+6
    return out

# ionization parameter from the SFR surface density and gas surface density
def U_Sigmag_Sigmasfr(Sigma_sfr, Sigma_g):
    out= (1.7e+14)*(Sigma_sfr/(Sigma_g*Sigma_g))
    return out

# ionization column density
def Ni(U,Z):
    ND=Nd(Z)
    tau_sd = (1e+23*U)/ND
    NN=ND*np.log((1+tau_sd)/(1+(tau_sd/np.exp(1.0))))
    return NN

# depth where A_V = 1
def Nd(Z):
    Ndsolar=1.7e21 #cm^-2
    return Ndsolar/Z

# ====================
# stuff related to the penetration of the Lyman Werner photons (see Ferrara+2019)
# ====================
def w_of_D(Z):
    w=1.0/(1.0+ 0.9*(Z)**0.5)
    return w

def chi_of_U(Z):
    chi=8.7e+4*Z
    return chi

def chi_prime(U, Z):
    ww=w_of_D(Z)
    chi=chi_of_U(U)
    chiprime=ww*chi
    return chiprime

def NL(U,Z):
    ND=Nd(Z)
    NN=ND*np.log(1+chi_prime(U,Z))
    return NN

# ====================
# HI column in the ionized layer
# ====================
def NHIyi(U, Z):
    Ns=1e+23*U
    tau_s  = 2.7e+5 * U
    out = (Ns/tau_s) * (1.0 - Ni(U,Z)/Nd(Z))
    return out

def NHIy0(U, Z, N0):
    Ns=1e+23*U
    y0=N0/Ns
    tau_sd = 59.0 * U * Z
    tau_s  = 2.7e+5 * U
    out = (Ns/tau_s) * np.log(tau_sd/np.abs((np.exp(tau_sd*y0)-tau_sd-1.0)))
    return out

# === CII flux for the density bounded case
# === dens bounded case =================
def fcii_DB(n, Z, U, column, TPDR, THII):
    g2_cii=4.
    g1_cii=2.
    E12_158um = 0.0079 #eV
    A21_158um = 2.4e-6 #s^-1
    if(n<=3300):
        # rates:  
        ladmbdaCIII  = lambdaCIIIe(THII)
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

# === ion bounded case, with N0 < NL ====
def fcii_IB_N0(n, Z, U, column, TPDR, THII):
    g2_cii=4.
    g1_cii=2.
    E12_158um = 0.0079 #eV
    A21_158um = 2.4e-6 #s^-1
    # ionized column density
    N_i= Ni(U,Z)
    
    if(n<=3300):
        # rates:  
        ladmbdaCIII  = lambdaCIIIe(THII)
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

# === ion bounded case, with NL < N0 ====
def fcii_IB_NL(n, Z, U, column, TPDR, THII):
    g2_cii=4.
    g1_cii=2.
    E12_158um = 0.0079 #eV
    A21_158um = 2.4e-6 #s^-1
    # ionized column density
    N_i= Ni(U,Z)
    # NL column density
    N_L=NL(U, Z)
    
    if(n<=3300):
        # rates:  
        ladmbdaCIII  = lambdaCIIIe(THII)
        lambdaCII    = lambdaCIIh(TPDR)
        lambdaCII4   = lambdaCIIe(THII)


        # cii from neutral
        fcii_neutral = n*carbon_abundance(Z)*lambdaCII*(N_L - N_i)
    
        # cii from ionized layer
        fcii_ionized_IB=n*carbon_abundance(Z)*lambdaCII4*NHIyi(U, Z)
    
        out= fcii_neutral + fcii_ionized_IB
        
    else:
        
        LTE_pop_levels_PDR = (g2_cii/g1_cii)*np.exp((-1.602e-12*E12_158um)/(1.38065e-16*TPDR))
        LTE_pop_levels_HII = (g2_cii/g1_cii)*np.exp((-1.602e-12*E12_158um)/(1.38065e-16*THII))
        fcii_neutral = LTE_pop_levels_PDR * carbon_abundance(Z) * A21_158um * (1.602e-12*E12_158um)* (N_L - N_i)
        fcii_ionized_IB = LTE_pop_levels_HII * carbon_abundance(Z) * A21_158um * (1.602e-12*E12_158um)* NHIyi(U, Z)
        out= fcii_neutral + fcii_ionized_IB
        
    return out


# ======================================
# ======================================
# equation for SigmaCII ================
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

def sigma_cii_IB_NL(logn, Z, k, Sigma_sfr):
    n= 10**logn
    SS_g=Sigmag_of_Sigmasfr(Sigma_sfr, k)
    UU=U_Sigmag_Sigmasfr(Sigma_sfr, SS_g)
    column_density=(SS_g*10**22.0)/7.5e+7
    ff=fcii_IB_NL(n, Z, UU, column_density, 100., 1e+4)
    SS_CII=ff*2.474e+9
    return SS_CII

# =========================================
def sigma_cii(logn, Z, k, Sigma_sfr):
    
    SS_g=Sigmag_of_Sigmasfr(Sigma_sfr, k)
    UU=U_Sigmag_Sigmasfr(Sigma_sfr, SS_g)
    column_density=(SS_g*10**22.0)/7.5e+7
    N_i=Ni(UU,Z)
    N_L=NL(UU,Z)
    
    if(column_density<N_i):
        out = sigma_cii_DB(logn, Z, k, Sigma_sfr)
    elif(column_density<N_L):
            out=sigma_cii_IB_N0(logn, Z, k, Sigma_sfr)
    else:
            out=sigma_cii_IB_NL(logn, Z, k, Sigma_sfr)
    return out

# ======================================
# ======================================
# equation for SigmaCIII ===============
def sigma_ciii(logn, Z, k, Sigma_sfr):
    n=10**logn
    SS_g=Sigmag_of_Sigmasfr(Sigma_sfr, k)
    UU=U_Sigmag_Sigmasfr(Sigma_sfr, SS_g)
    column=(SS_g*10**22.0)/7.5e+7
    ff=fciii(n, Z, UU, 1e+4)
    SS_CIII=ff*2.474e+9 # 
    return SS_CIII
    
# ==== deviation from DL ================
# =======================================
def Delta(logn, Z, k, Sigma_sfr):
    out = np.log10(sigma_cii(logn, Z, k, Sigma_sfr))-np.log10(delooze_fit_resolved(Sigma_sfr))
    return out

# ==== deviation from HC ================
# =======================================
def Delta_HC(logn, Z, k, Sigma_sfr):
    out = sigma_cii(logn, Z, k, Sigma_sfr)/herrera_fit_resolved(Sigma_sfr)
    return out


# ====== SIGMA OIII lines . ========
def collisionstrengths_88(T):
    coll=[0.5814, 0.5005, 0.4866, 0.55240, 0.5648, 0.6007,0.6116]
    temp=[100, 500, 1000, 5000, 10000, 20000, 30000]
    return np.interp(T, temp, coll)

def collisionstrengths_52(T):
    coll=[1.036, 1.032, 1.072, 1.210, 1.330, 1.451, 1.499, 1.5]
    temp=[100, 500, 1000, 5000, 10000, 20000, 30000, 1e+5]
    return np.interp(T, temp, coll)

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