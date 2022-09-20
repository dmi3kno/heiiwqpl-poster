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
