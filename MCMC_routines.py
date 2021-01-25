import numpy as np
import analytical_equations as eqs

class galaxy_template:

  def __init__(self
               ,Sigma_SFR  = 2.0
               ,Sigma_CII  = 3.0e+7
               ,Sigma_OIII = 7.0e+7
               #
               ,rel_err_Sigma_CII      = 0.2
               ,rel_err_Sigma_OIII     = 0.2
               ,rel_err_Delta         = 0.2
               ):

    self.Sigma_SFR       = Sigma_SFR  # Msun/yr/kpc^2
    self.Sigma_CII       = Sigma_CII  # Lsun/kpc^2
    self.Sigma_OIII      = Sigma_OIII # Lsun/kpc^2

    self.rel_err_Sigma_CII     = rel_err_Sigma_CII  # relative error on the SigmaCII
    self.rel_err_Sigma_OIII    = rel_err_Sigma_OIII # relative error on the SigmaOIII
    self.rel_err_Delta        = rel_err_Delta     # relative error on the Delta

  def data_for_MCMC(self):
      y    = np.array([self.Deltagalaxy()
                      ,self.Sigma_CII
                      ,self.Sigma_OIII
                      ])
      yerr = np.array([self.Deltagalaxy()*self.rel_err_Delta
                      ,self.Sigma_CII*self.rel_err_Sigma_CII
                      ,self.Sigma_OIII*self.rel_err_Sigma_OIII
                      ])

      par  = self.Sigma_SFR

      return y, yerr, par

  def Deltagalaxy(self):
      from analytical_equations import delooze_fit_resolved
      out = np.log10(self.Sigma_CII) -np.log10(delooze_fit_resolved(self.Sigma_SFR))

      return out

  def print_info(self):

      print("Galaxy input data")
      #
      print("  Sigma_SFR         = ",self.Sigma_SFR , "Msun/yr"   )
      print("  Sigma_CII         = ",self.Sigma_CII , "Lsun/kpc^2")
      print("  Sigma_OIII        = ",self.Sigma_OIII, "Lsun/kpc^2")
      #
      print("  delta Sigma_SFR   = ",100.0*self.rel_err_Delta ,"%")
      print("  delta Sigma_CIII  = ",100.0*self.rel_err_Sigma_CII ,"%")
      print("  delta Sigma_OIII  = ",100.0*self.rel_err_Sigma_OIII ,"%")

class MC_model:

    def __init__(self

            , lognMIN =  0.5
            , lognMAX =  3.5
            , logkMIN = -1.0
            , logkMAX =  2.5
            , logZMIN = -1.5
            , logZMAX =  0.0

            ):

       # priors for the model
       self.lognMIN = lognMIN   # minimum log density [cm-3]
       self.lognMAX = lognMAX   # maximum log density [cm-3]
       self.logkMIN = logkMIN   # minimum log k_s
       self.logkMAX = logkMAX   # maximum log k_s
       self.logZMIN = logZMIN   # minimum log metallicity [Zsun]
       self.logZMAX = logZMAX   # maximum log metallicity [Zsun]

       assert self.lognMIN < self.lognMAX
       assert self.logkMIN < self.logkMAX
       assert self.logZMIN < self.logZMAX


    def lnprior(self,theta):
       logn, logZ, logk = theta

       #flat priors on logn, logk, and log Z

       ok =        self.lognMIN < logn < self.lognMAX
       ok = ok and self.logkMIN < logk < self.logkMAX
       ok = ok and self.logZMIN < logZ < self.logZMAX

       if  ok:
           out = 0.0
       else:
           out = -np.inf
       
       return out

    def model(self,theta, ssfr):

       logn, logZ, logk  = theta
       delta_galaxy      = eqs.Delta(logn=logn    , Z=10**logZ, k=10**logk, Sigma_sfr=ssfr)
       sigma_cii_galaxy  = eqs.sigma_cii(logn=logn, Z=10**logZ, k=10**logk, Sigma_sfr=ssfr)
       sigma_oiii_galaxy = eqs.Sigma_OIII88(logn=logn, Z=10**logZ, k=10**logk, Sigma_sfr=ssfr)
       
       return delta_galaxy,sigma_cii_galaxy,sigma_oiii_galaxy


    def lnlike(self,theta, y, yerr, ssfr):
       mod              = np.asarray(self.model(theta=theta, ssfr=ssfr))
       inv_sigma2       = 1.0/(yerr**2)
       out              = -0.5*(np.sum((y-mod)**2*inv_sigma2))
       return out

    def lnprob(self,theta, y, yerr, ssfr):
        lp = self.lnprior(theta)

        if not np.isfinite(lp):
           out = -np.inf
        else:
           out = lp + self.lnlike(theta,y, yerr, ssfr)

        return out

    def print_info(self):
       print("Priors")
       print("{:10} {:10} {:10}".format(self.lognMIN,"< log(n/cm^-3) <",self.lognMAX))
       print("{:10} {:10} {:10}".format(self.logkMIN,"< log(k_s)     <",self.logkMAX))
       print("{:10} {:10} {:10}".format(self.logZMIN,"< log(Z/Z_sun) <",self.logZMAX))



