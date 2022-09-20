tar_target(fld_ppc_draws, {
  d_grd <- seq(1,280, by=10)
  p_grd <- c(0.05, 0.5, 0.95)
  N <- 500
  #fld_pqr_draws <- posterior::as_draws_array(fld_pqr_fit) 
  
  fld_pqr_draws %>% 
    posterior::subset_draws(draw=sample(posterior::ndraws(fld_pqr_draws), N)) %>%
    #posterior::rename_variables(
    #  alph = "Q(v)"#,
    #  #bt = "Q(w)"
    #) %>% 
    posterior::as_draws_df() %>%
    as_tibble() %>% 
    crossing(p_grd) %>% 
    mutate(ds=list(d_grd),
           ys=pmap(list(p_grd, alph, bt, eta, k, dlt, ds), qrqf),
           tbl=pmap(list(draw=.draw, ps=p_grd, ds=ds, ys=ys), dplyr::bind_cols)) %>% 
    select(tbl) %>% 
    unnest(tbl) 
  
})
