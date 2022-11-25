stantargets::tar_stan_mcmc(mod_mtx, stan_files = "stan/MtX.stan",
               data = birdsdata, init = initfun_mtx,# pars = params,
               parallel_chains = 4, iter_sampling = 15000-5000, iter_warmup = 5000, thin = 10,
               seed=42, max_treedepth = 15,
              variables = c("beta", "omega", "mu_size", "sd_size"))
