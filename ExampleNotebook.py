# %%
import numpy as np
import analytical_equations as eqs
from   MCMC_routines import MC_model,galaxy_template
import emcee

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

# %%
"""
# Set the ranges for the flat priors in the MCMC routines 
The ranges are bound between (lognMIN, lognMAX) (logZMIN, logZMAX), (logkMIN, logkMAX), with the following defaults
 - lognMIN = 0.5, lognMAX= 3.5
 - logZMIN = -1.5, logZMAX= 0
 - logkMIN = -1, logkMAX= 2.5

"""

mcr = MC_model(
  lognMIN = 0.5, lognMAX= 3.5
 ,logZMIN = -1.5, logZMAX= 0
 ,logkMIN = -1, logkMAX= 2.5
        )

# %%
# should you want to modify them, uncomment the following lines and set your preferred values
mcr.print_info()


# %%
"""
# Set up the MCMC details

For further details on the meaning of these parameters, 
please have a look at the emcee documentation at https://emcee.readthedocs.io/en/stable/
- <b> n_dim </b> i.e. the number of dimensions (in this case 3 because our model has three free parameters)
- <b> n_walkers </b> i.e. the the number of walkers
- <b> steps </b> i.e. the number of steps for each walker
- <b> burn_in </b> i.e. the number of initial steps that one may want to discard (the so-called "burn-in")
- <b> starting point </b>  i.e. the log$n_0$,log$Z_0$,log$k_0$ starting points around which the walkers are initialized
"""

# %%
n_dim               = 3 
n_walkers           = 10 
steps               = 200
burn_in             = 50
logn0, logZ0, logk0 = 2.0, -0.5, 0.3
starting_point      = [logn0, logZ0, logk0]
pos                 = [starting_point + 1e-5*np.random.randn(n_dim) for i in range(n_walkers)]

# %%
"""
# Run the MCMC
Further details on the possibile optimization of the MCMC algorithm and the outputs can be found in the documentation of emcee.
"""

# %%

y, yerr, par = galaxy_example.data_for_MCMC()

sampler      = emcee.EnsembleSampler(n_walkers, n_dim, mcr.lnprob, args=(y, yerr, par))
sampler.run_mcmc(pos, steps, progress=True)
tau          = sampler.get_autocorr_time(quiet=True)
flat_samples = sampler.get_chain(discard=burn_in, flat=True)

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



