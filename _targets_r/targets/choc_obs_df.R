tar_target(choc_obs_df, {
  set.seed(42)
  chocolate_nested %>% 
    filter(Survey=="(Adolescents, 2016)") %>% 
    #mutate(Ncons=100) %>% #Ncons*10
    mutate(obs=map2(Ncons, mtlg_a, qpd::rmetalog, bl=0))
})
