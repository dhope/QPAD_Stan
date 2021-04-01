data{
  // Data dimensions
  int<lower=1> n_sites;     // Number of surveys
  int<lower=1> n_ints;      // Number of censored detection intervals
  int<lower=0> y[n_sites,n_ints];  // Counts by survey and observation interval
  
// Covariates
// int<lower=0> n_bd;        // Number of detection fixed effects
// matrix[n_sites,n_bd] Xd;  // Explanatory variables for detection
  // Data
  real<lower=0> tau[n_ints];       // Endpoints of observation intervals
}

transformed data {
int<lower=0> ii[n_sites];         // An index identifying where in the vectorized data each survey begins (minus 1)
int<lower=0> yv[n_sites*n_ints];  // Vectorized data
int<lower=0> obsN[n_sites];       // Total count at each survey across all intervals

// Calculates each of the above
for (s in 1:n_sites) {
  ii[s] = (s-1)*n_ints;
  obsN[s] = 0;
  for (i in 1:n_ints) {
    yv[ii[s]+i] = y[s,i];
    obsN[s] = obsN[s] + y[s,i];
  }
}
}

parameters{
  real intcpt_d;     // Detection intercept
  // vector[n_bd] bd;   // Detection fixed effect estimates
  real<lower=0> sigma_det;        // Shape parameter for the Lognormal distribution
  real<lower=0,upper=1> gamma;    // Heterogeneity mixing parameter
}
transformed parameters {
  vector[n_sites] mu_det;                 // Survey-level mu from the LogNormal distribution (= -log(phi))
  vector<upper=0>[n_ints] log_p[n_sites]; // Interval-specific log(detection probability)
  vector[n_sites*n_ints] log_mu_ab;       // Interval-specific log(Expected count = lambda * Pr(detection in interval))

  // Calculating Survey-level rate of detection
  // mu_det = -Xd*bd;  // Note: the negative sign is in keeping with the other distributions, since E[T] \propto 1/rate
  for (s in 1:n_sites) {
    mu_det[s] =  intcpt_d;// mu_det[s] -- rd[id[s]]; // WE COULD ADD 0.5*(sigma_det)^2 TO MODEL THE TRUE MEAN
  }
  for (s in 1:n_sites) {
  log_p[s,1] = log(gamma * lognormal_cdf(tau[1],mu_det[s],sigma_det) + 1 - gamma);
  log_mu_ab[ii[s]+1] = log(obsN[s])+log_p[s,1];
  for (i in 2:n_ints) {
      log_p[s,i] = log(gamma) + log_diff_exp(lognormal_lccdf(tau[i-1]|mu_det[s],sigma_det),
                                              lognormal_lccdf(tau[i]|mu_det[s],sigma_det));
       log_mu_ab[ii[s]+i] = log(obsN[s])+log_p[s,i];
  }
  }
}
model{
  // Fixed effect priors
  // bd ~ normal(0,1.2); 
  intcpt_d ~ normal(-2.3,1.2);
  sigma_det ~ cauchy(0,1);
  gamma ~ beta(1,1);


  // Data model
  yv ~ poisson_log(log_mu_ab);
  }

generated quantities{
  // vector[n_sites*n_ints] log_lik;
  // real phi;
  matrix[n_sites,n_ints] phi_est_i;
  vector[n_sites] phi;
  for(s in 1:n_sites){
  phi_est_i[s,1] = exp(log_p[s,1]);
  for(i in 2:n_ints){
    phi_est_i[s,i] = exp(log_p[s,i]) + phi_est_i[s,i-1];
  }
  phi[s] = exp(-1 .* mean(phi_est_i[s,]));
}


}
