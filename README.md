Target Markdown
================

Target Markdown is a powerful R Markdown interface for reproducible
analysis pipelines, and the chapter at
<https://books.ropensci.org/targets/markdown.html> walks through it in
detail. This R Markdown report the example from the chapter. Try it out
in both interactive and non-interactive modes, either by running the
code chunks in different ways or setting the `tar_interactive` chunk
option.

# Packages

The example requires several R packages, and `targets` must be version
0.5.0.9000 or above.

``` r
remotes::install_github("dmi3kno/qpd")
install.packages(c("targets", "tarchetypes","extraDistr", "tidyverse", "posterior", "rstan", "fmcmc", "details"))
```

# Setup

If you are using old versions of `targets` (\<= 0.7.0) and/or `knitr`
(\<= 1.33), you will need to load the `targets` package in the R
Markdown document in order for Target Markdown code chunks to work.

``` r
library(targets)
library(magrittr)
tar_unscript()
```

# Globals

We first define some global options/functions common to all targets. The
function below plots a histogram of ozone concentrations, and our
histogram target will need it.

``` r
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("qpd", "extraDistr", "tidyverse", "posterior", "cmdstanr", "fmcmc"))
options(mc.cores = parallel::detectCores()-2)
options(clustermq.scheduler="multicode")
#> Establish _targets.R and _targets_r/globals/heiiwqplposter-globals.R.
```

# Targets

Prepare EFSA data

``` r
tar_target(chocolate_nested, {
  chocolate_raw <- readxl::read_excel("data/Foodex 2 L5 dashboard.xlsx", skip = 2) 
  chocolate <- chocolate_raw %>% 
    filter(`Survey start year`>2000) %>% 
    mutate(Survey=paste0("(", `Population Group`,", ",`Survey start year`,")")) %>% 
    select(21,10,14:19) %>% 
    rename("50th percentile"="Median",
           "Ncons"="Number of consumers")
  
  chocolate %>% 
    pivot_longer(-c(Survey,Ncons), names_to = "p", values_to = "q") %>% 
    mutate(p=parse_number(p)/100) %>% 
    group_by(Survey,Ncons,q) %>% 
    summarize(p=mean(p)) %>% 
    group_by(Survey, Ncons) %>% 
    summarise(p=list(p), q=list(q)) %>% 
    mutate(mtlg_a=map2(p, q, ~qpd::approx_max_metalog(q=.y, p_grid =.x, bl=0))) %>% 
    mutate(mtlg_v=map_lgl(mtlg_a, qpd::is_metalog_valid, bl=0))
})
```

Sample from latest adult survey metalog

``` r
tar_target(choc_obs_df, {
  set.seed(42)
  chocolate_nested %>% 
    filter(Survey=="(Adolescents, 2016)") %>% 
    #mutate(Ncons=100) %>% #Ncons*10
    mutate(obs=map2(Ncons, mtlg_a, qpd::rmetalog, bl=0))
})
#> Define target choc_obs_df from chunk code.
#> Establish _targets.R and _targets_r/targets/choc_obs_df.R.
```

Create a data frame of elicited data

``` r
tar_target(elicited_data, {
  tibble::tribble(~cat_idx, ~cat, ~CDFprob, ~count,
                    1L, "Morning coffee (0-4 g/day)", 0.25, 20,#18
                    1L, "Morning coffee (0-4 g/day)", 0.5, 24,#24
                    1L, "Morning coffee (0-4 g/day)", 0.75, 28,#30
                    2L, "After dinner (4-8 g/day)", 0.25, 49-24, #49
                    2L, "After dinner (4-8 g/day)", 0.50, 55-24, #55
                    2L, "After dinner (4-8 g/day)", 0.75, 60-24, #60
                    4L, "Sweet tooth (16+ g/day)", 0.25, 100-95,#92
                    4L, "Sweet tooth (16+ g/day)", 0.50, 100-85,#85
                    4L, "Sweet tooth (16+ g/day)", 0.75, 100-80,#80
    )
})
#> Define target elicited_data from chunk code.
#> Establish _targets.R and _targets_r/targets/elicited_data.R.
```

Fit a Dirichlet distribution to elicited data

``` r
tar_target(dir_df, {
  elicited_data %>% 
    mutate(prb=count/100) %>% 
    qpd::fit_dir(id_col = "cat_idx", spt_col = "CDFprob", prb_col="prb")
})
#> Define target dir_df from chunk code.
#> Establish _targets.R and _targets_r/targets/dir_df.R.
```

Prepare the data for QDirichlet-Metalog Stan model

``` r
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
#> Define target choc_data from chunk code.
#> Establish _targets.R and _targets_r/targets/choc_data.R.
```

Define initialization function for the QDir-Metalog model

``` r
initf_dirmetalog <- function() {
  deltas3 <- c(0.25, 0.25, 0.15)+runif(3, -0.02, 0.02)
  list(delta = c(deltas3, 1-sum(deltas3)))
}
#> Establish _targets.R and _targets_r/globals/f_initf_dirmetalog.R.
```

Compile and fit the model

``` r
# note the name of the target will be fit_dirmetalog_draws_mtqd_mod
stantargets::tar_stan_mcmc(name="fit_dirmetalog",
                       "stan/mtqd_mod.stan",
                       data=choc_data, init=initf_dirmetalog, 
         adapt_delta=0.8, seed=42, iter_sampling =5000, chains = 4)
#> Establish _targets.R and _targets_r/targets/model.R.
```

``` r
tar_target(draws_as_df, {
  m_as <- fit_dirmetalog_draws_mtqd_mod %>% 
    posterior::as_draws_array() %>% 
    posterior::subset_draws(variable="as", regex=TRUE)  %>% 
    posterior::as_draws_matrix() %>% unclass()
  
  map_df(seq_len(nrow(m_as)), 
         ~tibble::tibble(.draw=.x, as=list(m_as[.x,])))
})
#> Define target draws_as_df from chunk code.
#> Establish _targets.R and _targets_r/targets/draws_as_df.R.
```

# PQR

``` r
tar_target(carstop_data, {
  df_cs <- tibble::tibble(
    s=c(10, 10, 15, 15, 20, 20, 25, 25, 30, 30, 35, 35, 40, 40, 45, 45,
        50, 50, 55, 55, 60, 60, 65, 65, 70, 70, 75, 75, 80, 80),
    d= c(2, 2, 5, 6, 13, 11, 23, 21, 29, 36, 49, 50, 68, 71, 84, 81, 
         107, 99, 127, 132, 168, 122, 211, 195, 232, 176, 244, 263, 269, 236))
  df_cs
  
})
#> Define target carstop_data from chunk code.
#> Establish _targets.R and _targets_r/targets/carstop_data.R.
```

``` r
qrqf <- function(u, alph, bt, eta, k, dlt, x){
  S <- qpd::qfsld(u,eta,k,dlt)-qpd::qfsld(0.5, eta,k,dlt)
  alph+bt*sqrt(x)+sqrt(x)*S
}

irqf <- qpd::iqf(qrqf)

ll_rqf <- function(param, s, d){
  alph <- param[1]
  bt <- param[2]
  eta <- param[3]
  k <- param[4]
  dlt <- param[5]

  #if(v<0||v>1) return(-Inf)
  #if(w<0||w>1) return(-Inf)  
  if(eta<0) return(-Inf)
  if(k<0) return(-Inf)
  if(dlt<0||dlt>1) return(-Inf)

  # indirect priors
  #alph <- qpd::qfld(v, bt=1, k=10, a=1)
 # bt <- qpd::qfsld(w, bt=2, k=2, dlt=0.8, a=2)
  
  u <- irqf(s, alph, bt, eta, k, dlt, d, silent = FALSE)
# the scale adjustments need to be added to log-likelihood
  log_likelihood <- -log(sqrt(d))+qpd::dqfsld(u, eta, k, dlt, log=TRUE)
  #log_prior_alph <- dexp(alph, 1/5, log=TRUE) # quantile prior = 1/jacobian 
  log_prior_alph <- dlogitMyerson(bt, 0, 5, 11, alpha=0.1)
  #log_prior_bt <- dexp(bt, 1/5, log=TRUE) # quantile prior = 1/jacobian 
  log_prior_bt <- dlogitMyerson(bt, 2, 5, 12, alpha=0.1)
  log_prior_eta <- dexp(eta, 1/2, log=TRUE)
  log_prior_k <- dexp(k, 1/0.1, log=TRUE)
  log_prior_dlt <- dbeta(dlt, 2,1, log=TRUE)
  
  ll <- sum(log_likelihood)+log_prior_alph+
    log_prior_bt+
    log_prior_eta+log_prior_k+log_prior_dlt
  if(!is.finite(ll)) return(-Inf)
  ll
}
set.seed(42)
initials_rqf <- matrix(rep(c(0.5, 0.5, 2, 0.1, 0.75), each=4)+rnorm(4*5, 0, 0.01), nrow=4, byrow=FALSE, 
                       dimnames = list(NULL, c("alph", "bt", "eta", "k", "dlt")))
#> Establish _targets.R and _targets_r/globals/fld_pqr_globs.R.
```

``` r
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
#> Define target fld_pqr_fit from chunk code.
#> Establish _targets.R and _targets_r/targets/fld_pqr_fit.R.
```

``` r
tar_target(fld_pqr_draws, {
  fld_pqr_fit %>% 
    posterior::as_draws_array() %>% 
    #posterior::mutate_variables(
      #"Q(v)"= qpd::qfld(v, bt=1, k=10, a=1),
      #"Q(w)"= qpd::qfsld(w, bt=2, k=2, dlt=0.8, a=2)
    #) %>% 
    posterior::subset_draws(variable=c("alph", "bt", "eta", "k", "dlt"))
})
#> Define target fld_pqr_draws from chunk code.
#> Establish _targets.R and _targets_r/targets/fld_pqr_draws.R.
```

Posterior draws for 5-50-95 quantile

``` r
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
#> Define target fld_ppc_draws from chunk code.
#> Establish _targets.R and _targets_r/targets/fld_ppc_draws.R.
```

## Execute all targets

The following code will remove the old logs, re-build the outdated
targets and re-render the final paper. As a result a new directory
\_targets will be created on your computer and the hashed version of the
results (in compressed format) will be stored there.

If you inspect the scripts in the `_targets_r/` you will see a list of
statements like `tar_target(name=claims_data, {...})`. These targets
will be ran each individually, as the sub-tasks.

``` r
unlink("logs", recursive = TRUE)
targets::tar_make()
#> ✔ skip target fit_dirmetalog_file_mtqd_mod
#> ✔ skip target chocolate_nested
#> ✔ skip target carstop_data
#> ✔ skip target elicited_data
#> ✔ skip target choc_obs_df
#> ✔ skip target fld_pqr_fit
#> ✔ skip target dir_df
#> ✔ skip target fld_pqr_draws
#> ✔ skip target choc_data
#> ✔ skip target fld_ppc_draws
#> ✔ skip target fit_dirmetalog_data
#> ✔ skip target fit_dirmetalog_mcmc_mtqd_mod
#> ✔ skip target fit_dirmetalog_draws_mtqd_mod
#> ✔ skip target fit_dirmetalog_summary_mtqd_mod
#> ✔ skip target fit_dirmetalog_diagnostics_mtqd_mod
#> ✔ skip target draws_as_df
#> ✔ skip pipeline: 0.154 seconds
```

The targets dependency graph helps your readers understand the steps of
the pipeline at a high level. You can review the graph of task
dependencies, including the status of individual nodes, with

``` r
targets::tar_visnetwork()
```

The nodes get invalidated if any modification is made of the node or any
of its parents (ancestors). The out-of-date nodes will be rerun next
time you run tar_make().

At this point, you can go back and run {targets} chunks in interactive
mode without interfering with the code or data of the non-interactive
pipeline.

``` r
tar_progress_summary()
#> # A tibble: 1 × 6
#>   skipped started built errored canceled since        
#>     <int>   <int> <int>   <int>    <int> <chr>        
#> 1      16       0     0       0        0 0.044 seconds
```

The results can be retrieved by the name of the subtask. For example,
this will retrieve the latest fit object for the QDirichlet-Metalog
model from the article.

``` r
tar_read(fld_pqr_fit)|>lapply(head)
#> [[1]]
#>          alph       bt       eta          k       dlt
#> 2501 3.918549 4.469005 0.3222705 0.05596943 0.8759073
#> 2502 3.910127 4.503589 0.3079356 0.15611309 0.8876271
#> 2503 3.910127 4.503589 0.3079356 0.15611309 0.8876271
#> 2504 3.910127 4.503589 0.3079356 0.15611309 0.8876271
#> 2505 3.910127 4.503589 0.3079356 0.15611309 0.8876271
#> 2506 3.981110 4.470833 0.2236376 0.18665040 0.9007127
#> 
#> [[2]]
#>          alph       bt       eta          k       dlt
#> 2501 3.943055 4.470632 0.3031995 0.04678034 0.8867095
#> 2502 3.814230 4.531507 0.3049362 0.03238356 0.9180802
#> 2503 3.814230 4.531507 0.3049362 0.03238356 0.9180802
#> 2504 3.814230 4.531507 0.3049362 0.03238356 0.9180802
#> 2505 3.814230 4.531507 0.3049362 0.03238356 0.9180802
#> 2506 3.814230 4.531507 0.3049362 0.03238356 0.9180802
#> 
#> [[3]]
#>          alph       bt      eta         k       dlt
#> 2501 3.649209 4.622401 0.267893 0.0943123 0.7010758
#> 2502 3.649209 4.622401 0.267893 0.0943123 0.7010758
#> 2503 3.649209 4.622401 0.267893 0.0943123 0.7010758
#> 2504 3.649209 4.622401 0.267893 0.0943123 0.7010758
#> 2505 3.649209 4.622401 0.267893 0.0943123 0.7010758
#> 2506 3.649209 4.622401 0.267893 0.0943123 0.7010758
#> 
#> [[4]]
#>          alph       bt       eta           k       dlt
#> 2501 3.821236 4.436101 0.2194156 0.061723942 0.9179218
#> 2502 3.821236 4.436101 0.2194156 0.061723942 0.9179218
#> 2503 3.844867 4.429312 0.2269168 0.002758101 0.9015468
#> 2504 3.844867 4.429312 0.2269168 0.002758101 0.9015468
#> 2505 3.790381 4.469048 0.2439252 0.046617704 0.8711560
#> 2506 3.790381 4.469048 0.2439252 0.046617704 0.8711560
```

# Session Info

``` r
sessioninfo::session_info()%>%
  details::details(summary = 'Session Info')
```

<details closed>
<summary>
<span title="Click to Expand"> Session Info </span>
</summary>

``` r

─ Session info ───────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.2.1 (2022-06-23)
 os       Ubuntu 20.04.5 LTS
 system   x86_64, linux-gnu
 ui       X11
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Europe/Stockholm
 date     2022-09-20
 pandoc   2.18 @ /usr/lib/rstudio/bin/quarto/bin/tools/ (via rmarkdown)

─ Packages ───────────────────────────────────────────────────────────────────
 package     * version date (UTC) lib source
 backports     1.4.1   2021-12-13 [1] CRAN (R 4.2.0)
 base64url     1.4     2018-05-14 [1] CRAN (R 4.2.0)
 callr         3.7.1   2022-07-13 [1] CRAN (R 4.2.1)
 cli           3.3.0   2022-04-25 [1] CRAN (R 4.2.0)
 clipr         0.8.0   2022-02-22 [1] CRAN (R 4.2.0)
 codetools     0.2-18  2020-11-04 [4] CRAN (R 4.0.3)
 data.table    1.14.2  2021-09-27 [1] CRAN (R 4.2.0)
 desc          1.4.1   2022-03-06 [1] CRAN (R 4.2.0)
 details       0.3.0   2022-03-27 [1] CRAN (R 4.2.0)
 digest        0.6.29  2021-12-01 [1] CRAN (R 4.2.0)
 ellipsis      0.3.2   2021-04-29 [1] CRAN (R 4.2.0)
 evaluate      0.15    2022-02-18 [1] CRAN (R 4.2.0)
 fansi         1.0.3   2022-03-24 [1] CRAN (R 4.2.0)
 fastmap       1.1.0   2021-01-25 [1] CRAN (R 4.2.0)
 glue          1.6.2   2022-02-24 [1] CRAN (R 4.2.0)
 htmltools     0.5.2   2021-08-25 [1] CRAN (R 4.2.0)
 httr          1.4.3   2022-05-04 [1] CRAN (R 4.2.0)
 igraph        1.3.1   2022-04-20 [1] CRAN (R 4.2.0)
 knitr         1.39    2022-04-26 [1] CRAN (R 4.2.0)
 lifecycle     1.0.1   2021-09-24 [1] CRAN (R 4.2.0)
 magrittr    * 2.0.3   2022-03-30 [1] CRAN (R 4.2.0)
 pillar        1.8.1   2022-08-19 [1] CRAN (R 4.2.1)
 pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.2.0)
 png           0.1-7   2013-12-03 [1] CRAN (R 4.2.0)
 processx      3.7.0   2022-07-07 [1] CRAN (R 4.2.1)
 ps            1.7.1   2022-06-18 [1] CRAN (R 4.2.1)
 purrr         0.3.4   2020-04-17 [1] CRAN (R 4.2.0)
 R6            2.5.1   2021-08-19 [1] CRAN (R 4.2.0)
 rlang         1.0.5   2022-08-31 [1] CRAN (R 4.2.1)
 rmarkdown     2.14    2022-04-25 [1] CRAN (R 4.2.1)
 rprojroot     2.0.3   2022-04-02 [1] CRAN (R 4.2.0)
 rstudioapi    0.13    2020-11-12 [1] CRAN (R 4.2.0)
 sessioninfo   1.2.2   2021-12-06 [1] CRAN (R 4.2.0)
 stringi       1.7.8   2022-07-11 [1] CRAN (R 4.2.1)
 stringr       1.4.0   2019-02-10 [1] CRAN (R 4.2.0)
 targets     * 0.12.0  2022-04-19 [1] CRAN (R 4.2.0)
 tibble        3.1.8   2022-07-22 [1] CRAN (R 4.2.1)
 tidyselect    1.1.2   2022-02-21 [1] CRAN (R 4.2.0)
 utf8          1.2.2   2021-07-24 [1] CRAN (R 4.2.0)
 vctrs         0.4.1   2022-04-13 [1] CRAN (R 4.2.0)
 withr         2.5.0   2022-03-03 [1] CRAN (R 4.2.0)
 xfun          0.31    2022-05-10 [1] CRAN (R 4.2.0)
 xml2          1.3.3   2021-11-30 [1] CRAN (R 4.2.0)
 yaml          2.3.5   2022-02-21 [1] CRAN (R 4.2.0)

 [1] /home/dm0737pe/R/x86_64-pc-linux-gnu-library/4.2
 [2] /usr/local/lib/R/site-library
 [3] /usr/lib/R/site-library
 [4] /usr/lib/R/library

──────────────────────────────────────────────────────────────────────────────
```

</details>

<br>
