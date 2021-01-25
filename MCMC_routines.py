import numpy as np
import analytical_equations as eqs


lognMIN=0.5
lognMAX=3.5
logkMIN=-1
logkMAX=2.5
logZMIN=-1.5
logZMAX=0.


def lnprior(theta):
    logn, logZ, logk = theta

    #flat priors on logn, logk
    if lognMIN < logn < lognMAX and logkMIN < logk < logkMAX and logZMIN < logZ < logZMAX:
        return 0.0
    return -np.inf

def model(logn, logZ, logk, ssfr):

    k=10**logk
    Z=10**logZ
    Sigma_sfr= ssfr

    delta_galaxy = eqs.Delta(logn, Z, k, Sigma_sfr)
    sigma_cii_galaxy = eqs.sigma_cii(logn, Z, k, Sigma_sfr)
    sigma_oiii_galaxy = eqs.Sigma_OIII88(logn, Z, k, Sigma_sfr)

    return delta_galaxy,sigma_cii_galaxy,sigma_oiii_galaxy


def lnlike(theta, y, yerr, ssfr):
    logn, logZ, logk = theta
    mod = np.asarray(model(logn, logZ, logk, ssfr))
    inv_sigma2 = 1.0/(yerr**2)
    return -0.5*(np.sum((y-mod)**2*inv_sigma2))


def lnprob(theta, y, yerr, ssfr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta,y, yerr, ssfr)