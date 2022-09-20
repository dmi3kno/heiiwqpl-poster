tar_target(fld_pqr_fit, {
  .ub <- .Machine$double.xmax
  .lb <- -.Machine$double.xmax
  fit_rqf <- fmcmc::MCMC(ll_rqf,
                         initial = initials_rqf,
                         s=carstop_data$s, d=carstop_data$d,
                         nsteps = 5000,
                         burnin = 2500,
                         nchains = 4L,
                         kernel = kernel_ram(eps=1e-3,
                                             lb=c(.lb,.lb,0,0,0),
                                             ub=c(.ub,.ub,.ub,.ub,1)),
                         multicore = FALSE,
                         progress = TRUE)
  fit_rqf
})
