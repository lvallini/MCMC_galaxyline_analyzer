# General purpose and references

In this repository you can find the Python code that will allow you to derive the gas density, gas metallicity, and deviations from the Kennicutt-Schmidt relation of a galaxy with known star formation rate surface density (Sigma_SFR) and [CII] and [OIII] surface brigthness. 

Details on the equations, and the rationale behind the implemention of this method are provided in the following papers:
 
  - Vallini et al. 2021, submitted 
  - <a href="https://ui.adsabs.harvard.edu/abs/2020MNRAS.495L..22V/abstract">Vallini et al 2020</a> 
  - <a href="https://ui.adsabs.harvard.edu/abs/2019MNRAS.489....1F/abstract">Ferrara et al 2019</a> 

# The repository contains the following files:

- <b> analytical_equation.py </b> gathers all the analytical equations for the calculation of the [CII] and [OIII] surface brightnesses. Require Pyneb.
- <b> MCMC_routines.py </b> gathers all the routines related to the implementation of the MCMC algorithm brightnesses.
- <b> ExampleNotebook.ipynb </b>, a Jupyter notebook exemplifying how to use the analytical equations, run the MCMC model, and plot the results.

# Requirements
The modules require numpy, Pyneb, emcee, corner, matplotlib, scipy and other standard python libraries.
