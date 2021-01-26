
import numpy as np

def delooze_fit(logsfr):
    return (logsfr + 7.06)/1.0

def delooze_fit_resolved(Sigma_sfr):
    logSigma_cii=(np.log10(Sigma_sfr) +6.99)/0.93
    return 10**logSigma_cii

def delooze_delta(Sigma_sfr,Sigma_cii):
    return np.log10(Sigma_cii) - np.log10(delooze_fit_resolved(Sigma_sfr))



