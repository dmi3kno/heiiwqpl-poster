tar_target(draws_as_df, {
  m_as <- fit_dirmetalog_draws_mtqd_mod %>% 
    posterior::as_draws_array() %>% 
    posterior::subset_draws(variable="as", regex=TRUE)  %>% 
    posterior::as_draws_matrix() %>% unclass()
  
  map_df(seq_len(nrow(m_as)), 
         ~tibble::tibble(.draw=.x, as=list(m_as[.x,])))
})
