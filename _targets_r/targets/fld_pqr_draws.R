tar_target(fld_pqr_draws, {
  fld_pqr_fit %>% 
    posterior::as_draws_array() %>% 
    #posterior::mutate_variables(
      #"Q(v)"= qpd::qfld(v, bt=1, k=10, a=1),
      #"Q(w)"= qpd::qfsld(w, bt=2, k=2, dlt=0.8, a=2)
    #) %>% 
    posterior::subset_draws(variable=c("alph", "bt", "eta", "k", "dlt"))
})
