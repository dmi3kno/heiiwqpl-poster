options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("qpd", "extraDistr", "tidyverse", "posterior", "cmdstanr", "fmcmc"))
options(mc.cores = parallel::detectCores()-2)
options(clustermq.scheduler="multicode")
