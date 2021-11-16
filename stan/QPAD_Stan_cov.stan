functions{

  real get_remove(real dur, real phi){
    real neg_phi = -1 * (phi);
    return 1 - exp( dur * neg_phi );
    // this is the problem. Need to not include durations that were not recorded.
  }

  real get_dist(real dist, real tau){
    return  1 - exp(-1 *  pow(dist / tau,2) ); 
      }
      
   matrix get_PPsum(matrix mat_, int[] n_ints, vector var_, int type_, int n_sites){
     int max_ndist = max(n_ints);
  
    matrix[n_sites, max_ndist] CP;
    matrix[n_sites, max_ndist-1] CP_m1;
    matrix[n_sites, max_ndist] P;
    matrix[n_sites, max_ndist] PPsum;
    vector [n_sites] Psum;
    P = rep_matrix(0, n_sites, max_ndist);
    CP = rep_matrix(0, n_sites, max_ndist);
    for(i in 1:n_sites){
        for(j in 1:n_ints[i]){
          if(type_ == 1) {CP[i,j] = get_remove(mat_[i,j], var_[i]);}
          if(type_ == 0) {CP[i,j] = get_dist(mat_[i,j], var_[i] );
         }
      }
    }
    // print(dims(CP));
    // print(dims(CP_m1));
    // print(max_ndist);
    CP_m1 = block(CP, 1, 1, n_sites, (max_ndist-1));
    P = CP - append_col(rep_vector(0, n_sites), CP_m1);
    for(i in 1:n_sites){
          Psum[i] = 0;
          for(j in 1:n_ints[i]){
            Psum[i] += P[i,j];
          }
          PPsum[i,] = P[i,] / (Psum[i]);
        }
          
        
        //append_col(rep_matrix(0, n_sites,1), CP[,-n_ints[i] ]);
        // Psum[i] = rep_row_vector(0, n_ints[i]) * P[i,];
        // PPsum = P ./ (Psum);
        return(PPsum);
    }
        
        
 

}


data{
  // Data dimensions
  int<lower=1> n_sites_dur;     // Number of surveys
  int<lower=1> n_sites_dist;     // Number of surveys
  int<lower=1> n_ints_dur;      // Maximum Number of censored detection intervals - duration
  int<lower=1> n_ints_dist;      // Maximum Number of censored detection intervals - Distance
  int<lower=1> n_dur_parm;    // Number of parameters for duration
  int<lower=1> n_dist_parm;   // Number of parameters for distance
  // int<lower=1> n_dist;
  int y_dur[n_sites_dur,n_ints_dur];  // Counts by survey and observation interval
  int y_dist[n_sites_dist,n_ints_dist];  // Counts by survey and observation interval
  int<lower=0, upper=n_ints_dur> n_obs_int_dur_by_site[n_sites_dur];  // For variable length intervals, n
  int<lower=0, upper=n_ints_dist> n_obs_int_dist_by_site[n_sites_dist];  // For variable length intervals, n

// Covariates
// int<lower=0> n_bd;        // Number of detection fixed effects
// matrix[n_sites,n_bd] Xd;  // Explanatory variables for detection
  // Data
  matrix[n_sites_dur,n_ints_dur] dur;
  matrix[n_sites_dist,n_ints_dist] dist;
  // real<lower=0> tau[n_ints];       // Endpoints of observation intervals
  
  matrix[n_sites_dur, n_dur_parm] X_dur;
  matrix[n_sites_dist, n_dist_parm] X_dist;
}
transformed data{
  int dur_vec[3] = {3,5,10};
  real q_A[3] = {0.5, 1, -99};
}

parameters{
   vector[n_dur_parm] log_phi;
   vector[n_dist_parm] log_tau;
} 

transformed parameters{
  vector[n_sites_dur] phi;
  vector[n_sites_dist] tau;
  vector[n_sites_dur] log_lik_dur;
  vector[n_sites_dist] log_lik_dist;
  

  // Duration variables
  matrix[n_sites_dur, n_ints_dur] PP_dur;
  // Distance variables
  matrix[n_sites_dist, n_ints_dist] PP_dis;
  
  phi = exp(X_dur * log_phi);
  tau = exp(X_dist * log_tau);
  
  // matrix mat_, int[] n_ints, real var_, int type_, int n_sites
  PP_dur = get_PPsum(dur, n_obs_int_dur_by_site, phi, 1, n_sites_dur);
  PP_dis = get_PPsum(dist, n_obs_int_dist_by_site, tau, 0, n_sites_dist);
  for(i in 1:n_sites_dur){
      vector[n_obs_int_dur_by_site[i]] v;
      v = to_vector(sub_row(PP_dur, i, 1, n_obs_int_dur_by_site[i]));
      log_lik_dur[i] = multinomial_lpmf(y_dur[i, 1:n_obs_int_dur_by_site[i]] | v  );
        // y_dur[i,1:n_obs_int_dur_by_site[i] ] ~ multinomial(to_vector(PP_dur[i,1:n_obs_int_dur_by_site[i]]));
    }
  for(j in 1:n_sites_dist){
        vector[n_obs_int_dist_by_site[j]] k;
      k = to_vector(sub_row(PP_dis, j, 1, n_obs_int_dist_by_site[j]));
      log_lik_dist[j] = multinomial_lpmf(y_dist[j, 1:n_obs_int_dist_by_site[j]] | k  );
      }
}

model{
  
    log_phi ~ normal(0,0.5);//#beta(1,1);
    log_tau ~ normal(0,0.5);//#beta(1,1);
    target+=log_lik_dur;
    target+=log_lik_dist;
    

    
    
}
// generated quantities{
//   vector[3] availability_p;
//   vector[3] avg_detect;
//   vector[3] area_sampled;
//   vector[9] corr;
//   vector[9] off;
// 
//   for (i in 1:3){
//     availability_p[i] = 1 - exp(-1*dur_vec[i] * phi);
//     if(i<3)    {
//       avg_detect[i] = pow(tau,2) * (1-exp(-1* pow(q_A[i],2) / pow(tau,2) ) ) / pow(q_A[i],2);
//       area_sampled[i] = pow(q_A[i],2) * pi();
//     } else{
//       avg_detect[i] = 1;
//       area_sampled[i] = pow(tau,2) * pi();
//     }
//     }
// 
//     for(z in 1:3){
//       for(j in 1:3){
//           corr[ (z-1)*3 + j] = area_sampled[z]*availability_p[j]*avg_detect[z];
//         }
//     }
// 
//     off = log(corr);
// 
// 
// 
// 
// 
// 
// }
