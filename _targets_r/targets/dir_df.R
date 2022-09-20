tar_target(dir_df, {
  elicited_data %>% 
    mutate(prb=count/100) %>% 
    qpd::fit_dir(id_col = "cat_idx", spt_col = "CDFprob", prb_col="prb")
})
