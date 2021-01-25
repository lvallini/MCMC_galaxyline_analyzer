import numpy as np
import analytical_equations as eqs


class MC_model:


    def __init__(self

            , lognMIN=0.5
            , lognMAX=3.5
            , logkMIN=-1
            , logkMAX=2.5
            , logZMIN=-1.5
            , logZMAX=0.

            ):

       self.lognMIN = lognMIN
       self.lognMAX = lognMAX
       self.logkMIN = logkMIN
       self.logkMAX = logkMAX
       self.logZMIN = logZMIN
       self.logZMAX = logZMAX

    def lnprior(self,theta):
       logn, logZ, logk = theta

       #flat priors on logn, logk
       if self.lognMIN < logn < self.lognMAX and self.logkMIN < logk < self.logkMAX and self.logZMIN < logZ < self.logZMAX:
           out = 0.0
       else:
           out = -np.inf
       
       return out

    def model(self,logn, logZ, logk, ssfr):

       k=10**logk
       Z=10**logZ
       Sigma_sfr= ssfr

       delta_galaxy      = eqs.Delta(logn, Z, k, Sigma_sfr)
       sigma_cii_galaxy  = eqs.sigma_cii(logn, Z, k, Sigma_sfr)
       sigma_oiii_galaxy = eqs.Sigma_OIII88(logn, Z, k, Sigma_sfr)
       
       return delta_galaxy,sigma_cii_galaxy,sigma_oiii_galaxy


    def lnlike(self,theta, y, yerr, ssfr):
       logn, logZ, logk = theta
       mod = np.asarray(self.model(logn, logZ, logk, ssfr))
       inv_sigma2 = 1.0/(yerr**2)
       return -0.5*(np.sum((y-mod)**2*inv_sigma2))

    def lnprob(self,theta, y, yerr, ssfr):
        lp = self.lnprior(theta)
        if not np.isfinite(lp):
           return -np.inf
        return lp + self.lnlike(theta,y, yerr, ssfr)


    def print_info():
       print("Priors")
       print(self.lognMIN,"<log(n/cm^3)  <",self.lognMAX)
       print(self.logkMIN,"<log(k_s)     <",self.logkMAX)
       print(self.logZMIN,"<log(Z/Z_sun) <",self.logZMAX)



