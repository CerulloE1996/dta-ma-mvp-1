
functions {

          real corr(int num, int[] vec1, int[] vec2) {
            real pearson_corr =      (num*sum(to_vector(vec1) .* to_vector(vec2))   -  sum(to_vector(vec1))*sum(to_vector(vec2)) )  /  
                               sqrt( (num*sum(to_vector(vec1) .* to_vector(vec1))   -  sum(to_vector(vec1))^2) * (num*sum(to_vector(vec2) .*  to_vector(vec2))   - sum(to_vector(vec2))^2 ) ); 
            return(pearson_corr);
          }
real phi_logit_approx(real b) {
          real a;
          a = inv_logit( 1.702*b );
            return(a);
          }
  vector phi_logit_approx_vec(vector b) {
          int N = num_elements(b);
          vector[N] a;
          a = inv_logit( 1.702*b );
            return(a);
          }
  real inv_phi_logit_approx(real b) {
          real a;
          a = (1/1.702)*logit(b);    
            return(a);
          }


}

data {
          int num_binary_tests;
          int<lower=1> nt; // number of tests
          int<lower=1> Thr[nt]; // number of thresholds for each test
          int<lower=1> n_studies;	// total # of studies
          int<lower=1> ns[n_studies]; // number of individuals in each study 
          int<lower=0> y[max(ns),nt, n_studies]; // N individuals and nt tests, n_studies studies 
          int r[choose(nt,2), n_studies, 4];
          int numg;
          int total_n;
          int ns_cumsum[n_studies];
          int ind[n_studies];  
}


parameters {
              vector[nt] mu[2];  
              vector<lower=0>[2] sd[nt];
              vector[6] z1[n_studies];
              real<lower=0,upper=1> p[n_studies];   
              real<lower=0,upper=1> u_d[total_n,nt]; // nuisance that absorbs inequality constraints
              real<lower=0,upper=1> u_nd[total_n,nt]; // nuisance that absorbs inequality constraints
              cholesky_factor_corr[2] L_Omega_bs[nt];
              matrix[2,nt] z[n_studies];
}

transformed parameters {                         
                          matrix[2,nt] nu[n_studies];
                          vector[total_n] log_lik;
 
                      // between-study model
                for (s in 1:n_studies) 
                    for (t in 1:nt) 
                         nu[s, ,t] =  to_vector(mu[,t]) + diag_pre_multiply(sd[t,], L_Omega_bs[t,,]) * z[s, ,t] ;
      
                      // Parameters for likelihood function. Based on code upoaded by Ben Goodrich which uses the 
                      // GHK algorithm for generating TruncMVN. See:
                      // https://github.com/stan-dev/example-models/blob/master/misc/multivariate-probit/probit-multi-good.stan#L11
                      {
                        for (s in 1:n_studies) {
                              vector[nt] bound_d_bin;
                              vector[nt] bound_nd_bin;
                              bound_d_bin  = phi_logit_approx_vec(  -( to_vector(nu[s,1,1:nt]) ) );
                              bound_nd_bin = phi_logit_approx_vec(  -( to_vector(nu[s,2,1:nt]) ) );
                          for (n in 1:ns[s]) {
                              real lp1;
                              real lp0;
                              vector[nt] z_d;
                              vector[nt] z_nd;
                              vector[nt] y1d;
                              vector[nt] y1nd;
                            for (i in 1:nt) { // loop over each test - make sure binary tests are at the start 
                              real u_d1 =   u_d[n + ns_cumsum[s-1]*ind[s] ,i];
                              real u_nd1 = u_nd[n + ns_cumsum[s-1]*ind[s] ,i];
                              if (i <= num_binary_tests) { 
                                if (y[n,i,s] == 1) {
                                  z_d[i]   = inv_phi_logit_approx(bound_d_bin[i] + (1 - bound_d_bin[i])*u_d1);      
                                  z_nd[i]  = inv_phi_logit_approx(bound_nd_bin[i] + (1 - bound_nd_bin[i])*u_nd1);    
                                  y1d[i]   = log1m(bound_d_bin[i]);  // Jacobian adjustment
                                  y1nd[i]  = log1m(bound_nd_bin[i]); // Jacobian adjustment
                                }
                                if (y[n,i,s] == 0) { 
                                  z_d[i]   = inv_phi_logit_approx(bound_d_bin[i]*u_d1);
                                  z_nd[i]  = inv_phi_logit_approx(bound_nd_bin[i]*u_nd1); 
                                  y1d[i]   = log(bound_d_bin[i]);  // Jacobian adjustment
                                  y1nd[i]  = log(bound_nd_bin[i]); // Jacobian adjustment
                                }
                              }
                             // Wells
                              else { 
                                  if (y[n,i,s] == 2 || y[n,i,s] == 3) {  // moderate + high (2nd dichot)
                               //   if ( y[n,i,s] == 3) { // high (1st dichot)
                                  z_d[i]   = inv_phi_logit_approx(bound_d_bin[i] + (1 - bound_d_bin[i])*u_d1);      
                                  z_nd[i]  = inv_phi_logit_approx(bound_nd_bin[i] + (1 - bound_nd_bin[i])*u_nd1);    
                                  y1d[i]   = log1m(bound_d_bin[i]);  // Jacobian adjustment
                                  y1nd[i]  = log1m(bound_nd_bin[i]); // Jacobian adjustment
                                  }
                                if (y[n,i,s] == 1) {  // low
                             //   if (y[n,i,s] == 1 || y[n,i,s] == 2) { // low + moderate
                                  z_d[i]   = inv_phi_logit_approx(bound_d_bin[i]*u_d1);
                                  z_nd[i]  = inv_phi_logit_approx(bound_nd_bin[i]*u_nd1);
                                  y1d[i]   = log(bound_d_bin[i]);  // Jacobian adjustment
                                  y1nd[i]  = log(bound_nd_bin[i]); // Jacobian adjustment
                                }
                              }
                              // Jacobian adjustments imply z is truncated standard normal
                              // thus utility --- mu + L_Omega * z --- is truncated multivariate normal
                            }
                            lp1  = sum(y1d)  +   bernoulli_lpmf(1 | p[s]); 
                            lp0  = sum(y1nd) +   bernoulli_lpmf(0 | p[s]);
                           log_lik[n + ns_cumsum[s-1]*ind[s]] =  log_sum_exp(lp1 , lp0); 
                          }
                        }
                      }
}
model {
          for (t in 1:nt) 
            sd[t,] ~ normal(0, 0.5);
      
          for (t in 1:nt) 
            L_Omega_bs[t,,] ~ lkj_corr_cholesky(2);
           
          mu[2,1] ~ normal(- 1.7, 0.40); 
          mu[1,1] ~ normal(0.75, 0.40); 
          mu[, 2] ~ std_normal(); 
          
          for (s in 1:n_studies) {
           z1[s,] ~ std_normal(); 
           for (t in 1:nt) 
           z[s,1:2,t] ~ std_normal(); 
          }
       
          for (s in 1:n_studies) { // likelihood
            for (n in 1:ns[s]) 
               target +=  log_lik[ n + ns_cumsum[s-1]*ind[s] ];  
          }
  }

generated quantities {
    vector<lower=0,upper=1>[nt] Se; 
    vector<lower=0,upper=1>[nt] Sp; 
    vector<lower=0,upper=1>[nt] se[n_studies]; 
    vector<lower=0,upper=1>[nt] sp[n_studies]; 
    vector<lower=0,upper=1>[nt] fp[n_studies]; 
    matrix[n_studies,4] pr[choose(nt, 2)];  
    matrix[n_studies,4] e[choose(nt, 2)]; 
    matrix[n_studies,4] o[choose(nt, 2)]; 
    matrix[n_studies,4] ot[choose(nt, 2)]; 
    matrix[n_studies,4] dt[choose(nt, 2)];
    matrix[choose(nt, 2),n_studies] ec; 
    matrix[choose(nt, 2),n_studies] oc; 
    matrix[choose(nt, 2),n_studies] dc; 
    real Wells_DDimer_BTN_Se; 
    real Wells_DDimer_BTN_Sp; 
    real Wells_DDimer_BTP_Se; 
    real Wells_DDimer_BTP_Sp; 
    
      // Global Estimates
    Se[1]  =    phi_logit_approx(   mu[1,1]  );
    Sp[1]  =    phi_logit_approx( - mu[2,1]  );
    Se[2]  =    phi_logit_approx(   mu[1,2]  );
    Sp[2]  =    phi_logit_approx( - mu[2,2]  );
    Se[3]  =    phi_logit_approx(   mu[1,3]  );
    Sp[3]  =    phi_logit_approx( - mu[2,3]  );

  // estimate joint test accuracy estimates for Wells and D-Dimer
      Wells_DDimer_BTN_Se = Se[2]*Se[3] ;
      Wells_DDimer_BTN_Sp = 1- ( (1-Sp[2])*(1-Sp[3])  ); 
      Wells_DDimer_BTP_Se = 1- ( (1-Se[2])*(1-Se[3])  ); 
      Wells_DDimer_BTP_Sp = Sp[2]*Sp[3] ;

    {
      // Study-specific Estimates
      for (s in 1:n_studies)  {
           
        // study-specific accuracy estimates 
        se[s,1]   =      phi_logit_approx(  nu[s,1,1] );
        sp[s,1]   =      phi_logit_approx( -nu[s,2,1] );
        se[s,2]   =      phi_logit_approx(  nu[s,1,2] );
        sp[s,2]   =      phi_logit_approx( -nu[s,2,2] );
        se[s,3]   =      phi_logit_approx(  nu[s,1,3] );
        sp[s,3]   =      phi_logit_approx( -nu[s,2,3] );

        fp[s,] = 1-sp[s,];
        
       // ref vs d-dimer
        pr[1,s,1] =  p[s] * ( se[s,1] * se[s,2]         ) + (1-p[s]) * (  fp[s,1] * fp[s,2]          );  
        pr[1,s,2] =  p[s] * ( se[s,1] * (1-se[s,2])     ) + (1-p[s]) * (  fp[s,1]  * (1-fp[s,2])     );
        pr[1,s,3] =  p[s] * ( (1-se[s,1]) * se[s,2]     ) + (1-p[s]) * ( (1-fp[s,1]) * fp[s,2]       ); 
        pr[1,s,4] =  p[s] * ( (1-se[s,1]) * (1-se[s,2]) ) + (1-p[s]) * ( (1-fp[s,1]) * (1-fp[s,2])   );
       // ref vs wells
        pr[2,s,1] =  p[s] * ( se[s,1] * se[s,3]         ) + (1-p[s]) * (  fp[s,1] * fp[s,3]          );  
        pr[2,s,2] =  p[s] * ( se[s,1] * (1-se[s,3])     ) + (1-p[s]) * (  fp[s,1]  * (1-fp[s,3])     );
        pr[2,s,3] =  p[s] * ( (1-se[s,1]) * se[s,3]     ) + (1-p[s]) * ( (1-fp[s,1]) * fp[s,3]       ); 
        pr[2,s,4] =  p[s] * ( (1-se[s,1]) * (1-se[s,3]) ) + (1-p[s]) * ( (1-fp[s,1]) * (1-fp[s,3])   );
       // d-dimer vs wells
        pr[3,s,1] =  p[s] * ( se[s,2] * se[s,3]         ) + (1-p[s]) * (  fp[s,2] * fp[s,3]          );  
        pr[3,s,2] =  p[s] * ( se[s,2] * (1-se[s,3])     ) + (1-p[s]) * (  fp[s,2]  * (1-fp[s,3])     );
        pr[3,s,3] =  p[s] * ( (1-se[s,2]) * se[s,3]     ) + (1-p[s]) * ( (1-fp[s,2]) * fp[s,3]       ); 
        pr[3,s,4] =  p[s] * ( (1-se[s,2]) * (1-se[s,3]) ) + (1-p[s]) * ( (1-fp[s,2]) * (1-fp[s,3])   );
  
       
       // construct model-predicted 2-way table for each study
      for ( t in 1:choose(nt,2) ) {
        for (a in 1:4) { 
           ot[t,s,a] = r[t,s,a];
            e[t,s,a] = ns[s] * pr[t,s,a] ;                        // expected cell count (2x2 tables so 4)
            o[t,s,a] = ot[t,s,a] / ns[s] ;                         //# observed prob
           dt[t,s,a]    = ot[t,s,a] - e[t,s,a]; 
        }
          
        // Corr residuals (Qu et al, 1996)
        ec[t,s] = (pr[t,s,1] - (pr[t,s,1]+pr[t,s,2]) * (pr[t,s,1]+pr[t,s,3]) )/             // expected correlations
        sqrt((pr[t,s,1]+pr[t,s,2]) * (pr[t,s,1]+pr[t,s,3]) * (1-pr[t,s,1]-pr[t,s,2]) * (1-pr[t,s,1]-pr[t,s,3]));
        oc[t,s] = (o[t,s,1] - (o[t,s,1]+o[t,s,2]) * (o[t,s,1]+o[t,s,3])) /            // # observed correlations
          sqrt((o[t,s,1]+o[t,s,2]) * (o[t,s,1]+o[t,s,3]) * (1-o[t,s,1]-o[t,s,2]) * (1-o[t,s,1]-o[t,s,3]));
        dc[t,s] = oc[t,s] - ec[t,s];       //  # correlation residual
       }                                    
      }    
    }
  }

