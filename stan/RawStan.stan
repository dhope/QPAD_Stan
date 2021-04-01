data{
  // Data dimensions
  int<lower=1> n_sites;     // Number of surveys
  int<lower=1> n_ints;      // Number of censored detection intervals
  int<lower=0> y[n_sites,n_ints];  // Counts by survey and observation interval
  
// Covariates
// int<lower=0> n_bd;        // Number of detection fixed effects
// matrix[n_sites,n_bd] Xd;  // Explanatory variables for detection
  // Data
  matrix[n_sites,n_ints] tau;
  // real<lower=0> tau[n_ints];       // Endpoints of observation intervals
}
parameters{
   real<lower=0,upper=1> phi;
}

transformed parameters{
  matrix[n_sites, n_ints] CP;
  matrix[n_sites, n_ints] P;
  matrix[n_sites, n_ints] PPsum;
  vector[n_sites] RSP;
  real neg_phi;
  vector[n_sites] ll;
  neg_phi = -1 * phi;
  CP = 1 - exp( tau * neg_phi );
  P[,1] = CP[,1];
  for(i in 2:n_ints){
       P[,i] = CP[,i]-CP[,i-1]; 
  }
  for(j in 1:n_sites){
    RSP[j] = sum(P[j,]);
    for(i in 1:n_ints){
    PPsum[j,i] = P[j,i] /RSP[j];
    }
  }

  for (s in 1:n_sites){
    real tmpsum = 0;
    for(i in 1:n_ints){
      tmpsum += y[s,i] .* log(PPsum[s,i]) - lgamma(y[s,i] + 1);
    }
  ll[s] = lgamma(RSP[s] + 1) + tmpsum;
  }
}

model{
    phi ~ beta(1,1);
    target += ll;
}

