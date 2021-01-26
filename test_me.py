# %%
import numpy as np
from   MCMC_routines import MC_model,galaxy_template

# %%
"""
# Provide your input data:  
   - <b> Sigma_SFR  </b> (M$_{\odot}$/yr/kpc$^2$)
   - <b> Sigma_CII  </b> (L$_{\odot}$/kpc$^2$)
   - <b> Sigma_OIII </b> (L$_{\odot}$/kpc$^2$)

along with their relative errors
stadandard assumed error is 20%
"""

galaxy_example = galaxy_template(
                Sigma_SFR  = 2.0
               ,Sigma_CII  = 3.0e+7
               ,Sigma_OIII = 7.0e+7
              )
galaxy_example.print_info()

"""
change error for Sigma CII
"""

galaxy_example.set_relative_errors(rel_err_Sigma_CII = 0.1)
galaxy_example.print_info()

# %%
"""
# Set up the MCMC details

there are 3 kind of parameters:
    priors
    walker initial position
    MCMC hyperparameters
these can be set at initialization time of MC_model or at a later stage
default parameter are adopted if no initialization is provvided
"""

mcr = MC_model()

"""
Set the ranges for the flat priors in the MCMC routines 
The ranges are bound between (lognMIN, lognMAX) (logZMIN, logZMAX), (logkMIN, logkMAX), with the following defaults
"""

mcr.set_priors(
               lognMIN = 0.5, lognMAX= 3.5
              ,logZMIN = -1.5, logZMAX= 0
              ,logkMIN = -1, logkMAX= 2.5
              )

"""
<b> starting point </b>  i.e. the log$n_0$,log$Z_0$,log$k_0$ starting points around which the walkers are initialized
"""

mcr.set_walkers(
         logn0=2.0
        ,logZ0=-0.5
        ,logk0 = 0.3
        )

"""
set the MCMC parameters
For further details on the meaning of these parameters, 
please have a look at the emcee documentation at https://emcee.readthedocs.io/en/stable/
- <b> n_walkers </b> i.e. the the number of walkers
- <b> steps </b> i.e. the number of steps for each walker
- <b> burn_in </b> i.e. the number of initial steps that one may want to discard (the so-called "burn-in")
"""

mcr.set_mc_parameters(
         n_walkers = 10
        ,steps     = 200
        ,burn_in   = 50
        )

"""
print out your settings
"""
# %%
mcr.print_info()

"""
input the galaxy data you have constructed to the MCMC
"""

mcr.set_galaxy_data(galaxy_data = galaxy_example)

# %%
# %%


# %%
"""
# Run the MCMC
Further details on the possibile optimization of the MCMC algorithm and the outputs can be found in the documentation of emcee.
"""

# %%

flat_samples = mcr.run_model(verbose=True)

# %%
"""
# Plot the result
For more information and details please refer to: https://corner.readthedocs.io/en/latest/pages/quickstart.html
"""

# %%
import corner
fig = corner.corner(flat_samples, labels=["log(n/cm$^{-3}$)", "log(Z/Z$_{\odot}$)", "log($\kappa_s$)"])

fig.savefig("test.png")

# %%



