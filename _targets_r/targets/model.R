# note the name of the target will be fit_dirmetalog_draws_mtqd_mod
stantargets::tar_stan_mcmc(name="fit_dirmetalog",
                       "stan/mtqd_mod.stan",
                       data=choc_data, init=initf_dirmetalog, 
         adapt_delta=0.8, seed=42, iter_sampling =5000, chains = 4)
