// begin file stan/mtqd_mod.stan
functions{ // begin the functions block

//  //given x and a grid, returns a linearly interpolated y
vector vlookup(vector x, vector xs, vector ys){
  int N =rows(x); // number of ordered!! datapoints
  int M =rows(xs); // number of ordered!! grid cells
  real t;
  real w;
  vector[N] res;
  int i = 1;
  for (j in 2:M){
    while (i<=N && x[i]<=xs[j] && x[i]>=xs[j-1]) {
     res[i] = ys[j-1];
     i+=1;
    }
  }
  return res;
}

vector fit_logmetalog3(vector ps, vector qntls, real bl) {
  /* fit metalog to the qp-pairs. Number of metalog terms is automatically
  matched to the number of qp-pairs.
  @ ps cumulative probabilities
  @ args quantiles quantiles (values of x) corresponding to cumulative probabilities p
  @ returns returns of vector of metalog parameters ("a"-coefficients)
  */
    int n = rows(ps);
    int m = rows(qntls);
    matrix[n, n] Y;
    vector[n] log_odds = logit(ps);
    vector[n] pmhalf =  ps - 0.5;
    vector[m] z = log(qntls-rep_vector(bl, n));
    int odd = 1;
    if(n != 3 || m != 3) reject("Fit_metalog3() is a 3-term metalog function!");
    Y[ , 1] = rep_vector(1, n);
    Y[ , 2] = log_odds;
    Y[ , 3] = pmhalf .* log_odds;
    return Y \ z;
  }

real logmetalog_s_qf(real p, vector a, real bl){
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

vector logmetalog_reparameterize(vector theta, int[ ] idx, vector qntls, real bl){
  /* Change parameterization from the simplex parameterization to the Taylor series parameterization ("a"-coefficients)
  @ args theta p-simplex vector of probabilities of length N
  @ args idx integer array of indices correponding to p-simplex (of length N). This is data argument.
  @ args qntls values of x correponding to quantiles determining the metalog (length N-1). This is data argument.
  @ return vector of metalog parameters ("a"-coefficients)
  */
  int N = rows(theta);
  int M = rows(qntls);
  vector[N] cumtheta;
  vector[M] ps;
  if (N != M+1) reject("Size of quantiles should be 1 less than a size of the simplex!");
  if (size(idx) != N) reject("Simplex size did not match the index!");
  cumtheta = cumulative_sum(theta[idx]);
  ps = cumtheta[1:M];
  return fit_logmetalog3(ps, qntls, bl);
}


real metalog_s_qdf(real p, vector a){
  /* Quantile density function of metalog distribution
  @ args p real value of cumulative probability
  @ args a vector of metalog parameters ("a"-coefficients)
  @ return returns real value of variable x correponding to cumulative probability p
  */
    int n = rows(a);
    real res=0.0;
    real pt1mp = p*(1-p);
    real logitp = logit(p);
    real pmhalf = p-0.5;
    int odd = 1;
    if (n == 0) reject("a cannot be of size zero");
    if(p>0){
      if (n > 1) res += a[2]/pt1mp;
      if (n > 2) res += a[3]*(pmhalf/pt1mp+logitp);
      if (n > 3) res += a[4];
      if (n > 4) {
        for (m in 5:n) {
          if (odd) {
            res += a[m]*(m-1)/2.0*pow(pmhalf, (m-3)/2.0);
          } else res += a[m]*(pow(pmhalf, m/2.0-1)/pt1mp+(m/2.0-1)*pow(pmhalf, m/2.0-2)*logitp);
         odd = odd == 0;
         }
      }
    }
    return res;
 }


vector logmetalog3_v_qf(vector p, vector a, real bl){
  int N = rows(p);
  vector[N] res;
  for (i in 1:N) {
    if(metalog_s_qdf(p[i], a) < 0) reject("Invalid metalog, rejecting!");
    res[i] = logmetalog_s_qf(p[i],a, bl);}
  return res;
}

vector logmetalog_iqf_algebra_system(vector u0, vector a, data real[] x_r, data int[] x_i){
  return [x_r[1] - logmetalog_s_qf(u0[1], a, x_r[2])]';
}

real approx_cdf_algebra(data real x, real u_guess, vector a, data real bl, data real rel_tol, data real f_tol, data real max_steps){
  return algebra_solver(logmetalog_iqf_algebra_system, [u_guess]', a, {x, bl}, {0},  rel_tol, f_tol, max_steps)[1];
}

real logmetalog_s_ldqf_lpdf(real p, vector a, data real bl){
  // has _lpdf suffix so that Stan wouldn't complain about non-PDF likelihood
  real res;
  //rel_tol = 1e-10, f_tol = 1e-6, and max_steps = 1e3
  res = inv(metalog_s_qdf(p,a)*logmetalog_s_qf(p, a, 0.0));
  if (res<=0) reject ("LRQDF error: Quantile density returned", res, ". Rejecting!");
  return log(res);
  }

} // finish the functions block

data {
  int<lower=0> N; // size of the data
  vector[N] x; // data on variable level
  int<lower=0> M; // size of probability grid
  vector[M] ys_grd; // probability grid
  vector[4] a; // dirichlet parameter vector
  int idx[4]; // index into the dirichlet bins
  vector[3] qntls; // vector of quantile values
  real bl; // lower bound. Defaults to zero
  real  rel_tol; // for algebra solver
  real f_tol; // for algebra solver
  real max_steps;// for algebra solver
}

transformed data{
  vector[N] x_srt = sort_asc(x);
}

parameters {
  simplex[4] delta; // dirichlet sample - a simplex
}

transformed parameters{
  // go grom indexed-theta-quantile parameterization to metalog a-coeff parameterization
  vector[3] as = logmetalog_reparameterize(delta, idx, qntls, bl);
}

model {
  vector[N] u;
  // create grid of xs given the parameter
  vector[M] xs_grd = logmetalog3_v_qf(ys_grd, as, bl);
  vector[N] u_guess = vlookup(x_srt, xs_grd, ys_grd);

  //Grids are done. Sampling!
  target += dirichlet_lpdf(delta | a);
  for (i in 1:N){
   u[i] = approx_cdf_algebra(x_srt[i], u_guess[i], as, bl, rel_tol, f_tol, max_steps);
   target += logmetalog_s_ldqf_lpdf(u[i] | as, bl);
  }
}
// end file stan/mtqd_mod.stan
