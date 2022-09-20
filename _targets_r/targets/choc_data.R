tar_target(choc_data, {
  q_star <- c(4,8,16)
  prior_idx <- dir_df$cat_idx
  prior_dir_a <- dir_df$a
  choc_obs <- choc_obs_df %>% 
      unnest(cols="obs") %>% pull()
  # create probability grid to be used for interpolation
  ys_grd <- qpd::make_pgrid(1000,2)
  
  list(N=length(choc_obs), x=choc_obs, M=length(ys_grd), ys_grd=ys_grd,
       a=prior_dir_a, idx=prior_idx,
       qntls=q_star, bl=0.0,  rel_tol=1e-15, f_tol=1e-10, max_steps=1e2)
})
