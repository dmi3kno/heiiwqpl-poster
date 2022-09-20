// begin stan/inv_sampl_mod.stan
functions{

real logmetalog_qf_s_cdf(real p, vector a, real bl){
  // code by Dmytro Perepolkin <dperepolkin@gmail.com>
  /* Quantile function of metalog distribution
  @ args p real value of cumulative probability
  @ args a vector of metalog parameters ("a"-coefficients)
  @ return returns real value of variable x correponding to cumulative probability p
  */
  real res=bl;
  int n = rows(a);
  real logitp = logit(p);
  real pmhalf = p-0.5;
  int odd = 1;
  if(p>0){
    res += a[1];
    if (n == 0) reject("a cannot be of size zero");
    if (n > 1) res += a[2]*logitp;
    if (n > 2) res += a[3]*pmhalf*logitp;
    if (n > 3) res += a[4]*pmhalf;
    if (n > 4) {
      for (m in 5:n) {
        if (odd) {
          res += a[m]*pow(pmhalf, (m-1)/2.0);
        } else res += a[m]*pow(pmhalf, m/2.0-1)*logitp;
        odd = odd == 0;
      }
    }
    res=exp(res);
  }
  return res;
}

} // end of functions block
data {
  int<lower=0> N;
  real<lower=0> y[N];
  int<lower=2> M; // number of terms for metalog
  vector[M] prior_as;
  real prior_bl;
}
parameters {
  real<lower=0, upper=1> u;
}
transformed parameters{
  real lambda=logmetalog_qf_s_cdf(u, prior_as, prior_bl);
}
model {
  u ~ uniform(0,1);
  for (n in 1:N)
    y[n] ~ exponential(lambda);
}
// end stan/inv_sampl_mod.stan
