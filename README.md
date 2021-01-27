# General purpose and references

In this <a href='https://github.com/lvallini/MCMC_galaxyline_analyzer'>repository</a> you can find the Python code that will allow you to derive the gas density, gas metallicity, and deviations from the Kennicutt-Schmidt relation of a galaxy with known star formation rate surface density (Sigma_SFR) and [CII] and [OIII] surface brigthness. 

Details on the equations, and the rationale behind the implemention of this method are provided in the following papers:
 
  - Vallini et al. 2021, submitted 
  - <a href="https://ui.adsabs.harvard.edu/abs/2020MNRAS.495L..22V/abstract">Vallini et al. 2020</a> 
  - <a href="https://ui.adsabs.harvard.edu/abs/2019MNRAS.489....1F/abstract">Ferrara et al. 2019</a> 

# The repository contains the following files:

- <b> cooling_rates.py, emission_models.py, empirical.py, ion_structure.py </b> gather all the analytical equations (mostly from Ferrara et al. 2019) for the calculation of the [CII] and [OIII] surface brightnesses.

- <b> MCMC_routines.py </b> gathers all the routines related to the implementation of the Markov Chain Monte Carlo algorithm.

- <a href='https://github.com/lvallini/MCMC_galaxyline_analyzer/blob/main/ExampleNotebook-v2.ipynb'> ExampleNotebook-v2.ipynb</a>, a Jupyter Notebook exemplifying how to set your input data, run the MCMC model, and plot the resulting likelihood distribution for the gas density, gas metallicity, and deviation from the Kennicutt-Schmidt relation.

# Requirements
The modules **require** numpy, scipy, matplotlb, <a href="https://github.com/Morisset/PyNeb_devel">Pyneb</a>, 
<a href='https://emcee.readthedocs.io/en/stable'>emcee</a>,  <a href="https://corner.readthedocs.io/en/latest/index.html">corner</a>.

