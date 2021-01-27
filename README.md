# General purpose and references

<a href='https://github.com/lvallini/MCMC_galaxyline_analyzer'>In this repository</a> you can find the Python code that will allow you to derive the gas density, gas metallicity, and deviations from the Kennicutt-Schmidt relation of a galaxy with known star formation rate surface density (Sigma_SFR) and [CII] and [OIII] surface brigthness. 

Details on the rationale behind the model implemention and on the equations are discussed in the following papers:
 
  - Vallini et al. 2021, submitted 
  - <a href="https://ui.adsabs.harvard.edu/abs/2020MNRAS.495L..22V/abstract">Vallini et al. 2020</a> 
  - <a href="https://ui.adsabs.harvard.edu/abs/2019MNRAS.489....1F/abstract">Ferrara et al. 2019</a> 

# The repository contains the following files:

- <a href='https://github.com/lvallini/MCMC_galaxyline_analyzer/blob/main/ExampleNotebook-v2.ipynb'> ExampleNotebook-v2.ipynb</a>, a Jupyter Notebook exemplifying how to set your input data, run the MCMC model, and plot the resulting likelihood distribution for the gas density, gas metallicity, and deviation from the Kennicutt-Schmidt relation.

- <b> cooling_rates.py, emission_models.py, empirical.py, ion_structure.py </b> gather all the analytical equations (mostly from Ferrara et al. 2019) for the calculation of the [CII] and [OIII] surface brightnesses.

- <b> MCMC_routines.py </b> gathers all the routines related to the implementation of the Markov Chain Monte Carlo algorithm.

You can clone (or download) the entire repository and try the model with input data and priors of your choice.

# Requirements
The modules **require** numpy, scipy, matplotlib, <a href="https://github.com/Morisset/PyNeb_devel">Pyneb</a>, 
<a href='https://emcee.readthedocs.io/en/stable'>emcee</a>,  <a href="https://corner.readthedocs.io/en/latest/index.html">corner</a>.

## Acknowledging this code in Scientific Publications

<div class="row codice">
<pre><code><span>@ARTICLE{Ferrara:2019,
       author = <span>{</span> {Ferrara}, A. and {Vallini}, L. and {Pallottini}, A. and {Gallerani}, S. and {Carniani}, S.
                 and {Kohandel}, M. and {Decataldo}, D. and {Behrens}, C.},
        title = "{A physical model for [C II] line emission from galaxies}",
      journal = {\mnras},
         year = 2019,
        month = oct,
       volume = {489},
       number = {1},
        pages = {1-12},
          doi = {10.1093/mnras/stz2031},
archivePrefix = {arXiv},
       eprint = {1908.07536},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2019MNRAS.489....1F},
   }</span></code>
</pre>
</div>


<div class="row codice">
<pre><code><span>@ARTICLE{Vallini:2020,
       author = <span>{</span> {Vallini}, L. and {Ferrara}, A. and {Pallottini}, A. and {Carniani}, S. and {Gallerani}, S.},
        title = "{Star formation law in the epoch of reionization from [C II] and C III] lines}",
      journal = {\mnras},
     keywords = {photodissociation region (PDR), galaxies: high-redshift, galaxies: ISM, Astrophysics - Astrophysics of Galaxies},
         year = 2020,
        month = jun,
       volume = {495},
       number = {1},
        pages = {L22-L26},
          doi = {10.1093/mnrasl/slaa047},
archivePrefix = {arXiv},
       eprint = {2003.06443},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2020MNRAS.495L..22V},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}</span>
</code></pre>
</div>

## Funding
This work has received funding by the ERC
<img src="https://p7.hiclipart.com/preview/510/516/838/european-research-council-french-institute-for-research-in-computer-science-and-automation-university-of-paris-saclay-dan-zhang.jpg" width="400" height="400" alt="ERClogo">
