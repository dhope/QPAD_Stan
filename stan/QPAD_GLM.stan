functions{

  matrix get_remove(matrix dur, real phi){
    real neg_phi = -1 * (phi);
    return 1 - exp( dur * neg_phi );
  }

  matrix get_dist(matrix dist, real tau, int n_sites, int n_dist){
    matrix[n_sites, n_dist] dout;
    for(i in 1:n_sites){
      for(j in 1:n_dist){
        dout[i,j] =  1 - exp(-1 *  (dist[i,j] / tau) ^2); // 1 - exp(-(r/tau)^2)
      }
    }
  return dout;
  }

  matrix get_P(matrix CP,int n_sites, int n_ints){
    matrix[n_sites, n_ints] P;
     P[,1] = CP[,1];
    for(i in 2:n_ints){
       P[,i] = CP[,i]-CP[,i-1];
    }
    return P;
  }
  vector get_RSP(matrix P,int n_sites, int n_ints){
      vector[n_sites] RSP;
       for(j in 1:n_sites){
      RSP[j] = sum(P[j,]);
       }
       return RSP;
  }


  matrix get_PP(matrix P,vector RSP,int n_sites, int n_ints){
    matrix[n_sites, n_ints] PPsum;
    for(j in 1:n_sites){
      for(i in 1:n_ints){
        PPsum[j,i] = P[j,i] /RSP[j];
    }
    }
  return PPsum;
  }

  vector get_loglik(int[,] y, matrix PPsum, vector RSP, int n_sites, int n_ints){
    vector[n_sites] ll;
      for (s in 1:n_sites){
    real tmpsum = 0;
    for(i in 1:n_ints){
      tmpsum += y[s,i] .* log(PPsum[s,i]) - lgamma(y[s,i] + 1);
    }
  ll[s] = lgamma(RSP[s] + 1) + tmpsum;
    }
    return ll;
  }

  real calc_off(real tau, real phi, real dist, real time){
    real p;
    real q;
    real A;
      if(dist > 1e5){
          A = pi()*tau^2;
          q = 1;
      } else{
        A = pi()*dist^2;
        q = tau^2 * (1 - exp(-dist^2/tau^2))/dist^2;
      }
      p = 1 - exp(-time * phi);

      return log(A*p*q);

  }

}


data{
  // Data dimensions
  int<lower=1> n_sites;     // Number of surveys
  int<lower=1> n_ints;      // Number of censored detection intervals
  int<lower=1> n_dist;      // Number of distance rings
  int<lower=1> n_pred;      // Number of covariates
  int<lower=1> n_dat;       // number of data points for glm
  int<lower=0> y_dur[n_sites,n_ints];  // Counts by survey and observation interval
  int<lower=0> y_dist[n_sites,n_ints];  // Counts by survey and observation interval

// Covariates
// int<lower=0> n_bd;        // Number of detection fixed effects
// matrix[n_sites,n_bd] Xd;  // Explanatory variables for detection
  // Data
  matrix[n_sites,n_ints] dur;
  matrix[n_sites,n_ints] dist;
  // real<lower=0> tau[n_ints];       // Endpoints of observation intervals

  int y_count[n_dat];  // Counts for glm
  real time[n_dat]; // Time intervals for each points
  real glm_dist[n_dat]; // Max bands for glm points
  matrix[n_dat, n_pred] Xd; // Glm covariates


}
parameters{
   real<lower=0,upper=1> phi;
   real<lower=0,upper=1> tau;
   vector[n_pred] ba;
}

transformed parameters{

  // Duration variables
  matrix[n_sites, n_ints] CP_dur;
  matrix[n_sites, n_ints] PP_sum_dur;
  matrix[n_sites, n_ints] P_dur;
  vector[n_sites] RSP_dur;
  vector[n_sites] ll_dur;
  // Distance variables
  matrix[n_sites, n_ints] CP_dist;
  matrix[n_sites, n_ints] PP_sum_dis;
  matrix[n_sites, n_dist] P_dist;
  vector[n_sites] RSP_dist;
  vector[n_sites] ll_dist;
  
  vector[n_dat] off;
  vector[n_dat] lambda;

  CP_dur = get_remove(dur, phi);
  P_dur = get_P(CP_dur,n_sites, n_ints);
  RSP_dur = get_RSP(P_dur, n_sites,  n_ints);
  PP_sum_dur = get_PP(P_dur,RSP_dur, n_sites, n_ints);
  ll_dur = get_loglik(y_dur, PP_sum_dur, RSP_dur, n_sites, n_ints);

  CP_dist = get_dist(dist, tau, n_sites, n_dist);
  P_dist = get_P(CP_dist,n_sites, n_dist);
  RSP_dist = get_RSP(P_dist, n_sites,  n_dist);
  PP_sum_dis = get_PP(P_dist, RSP_dist, n_sites, n_dist);
  ll_dist = get_loglik(y_dist, PP_sum_dis, RSP_dist, n_sites, n_dist);

  for( i in 1:n_dat){
    off[i] = calc_off(tau, phi, glm_dist[i], time[i]);

  }
  lambda = Xd * ba + off;

}

model{
    phi ~ beta(1,1);
    tau ~ beta(1,1);
    ba ~ normal(0,1);
    y_count ~ poisson_log(lambda);
    target += ll_dur;
    target += ll_dist;
}

