//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// NOTE: this is for two biomarkers only!! 
functions {
  matrix to_triangular(vector x, int K) {
    // could check rows(y) = K * (K + 1) / 2
    matrix[K, K] y;
    int pos = 1;
    for (i in 1: K) {
      for (j in 1:i) {
        y[i, j] = x[pos];
        pos += 1;
      }
      for (j in (i+1):K) {
        y[i, j] = 0;
      }
    }
    return y;
  }
  
 matrix to_lower_tri(vector x, int K){
  matrix[K,K] A = rep_matrix(0, K, K); 
   int pos = 1 ;
   for(i in 1: K){
     for(j in 1:K){
       if(j<i){
         A[i,j] = x[pos];
         pos +=1;
       } else if (j==i){
         A[i,j] = 1;
         
       }
     }
   }
   return(A);
   }
}

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> N; //total number of observations 
  int<lower=1> I; // number of individuals 
  int<lower=1> ids[N];
  int<lower=1>Q; //no. of biomarkers 
 // int<lower=1> P; //no of mean basis functions
  int<lower=1> L;
  int<lower=1> L_S;
  matrix[N,Q] x;
  vector[N] time;
  vector[N] y;
  matrix[N, L] b_spline;
  matrix[N, L_S] S_b_spline;
//  matrix[Q,Q] Sigma0;
// input these as data for now; we need to estimate these later
  int<lower = 0, upper = 1> simulate; 
}

// The parameters accepted by the model.
parameters {
 // for the longitudinal submodel 
//  matrix<lower=0>[Q,L_S] theta_sigma;
  vector[L_S] S_beta[Q];
  matrix[L_S,Q] Theta_raw[I];  // prior correlation (for the covariance matrix of Bijs)
  //matrix[Q, L_S] Theta[I];
  vector<lower=0>[L] S_tau[Q];
  cholesky_factor_corr[Q] Omega;  // prior correlation (for the covariance matrix of Bijs)
  vector<lower=0>[Q] tau;
  cholesky_factor_corr[L] S_Omega[Q];  // prior correlation (for the covariance matrix of Bijs)
  vector<lower=0>[L_S] theta_tau[Q]; 
  cholesky_factor_corr[L_S] theta_Omega[Q];  // prior correlation (for the covariance matrix of Bijs)
  vector[L] beta[Q];
  matrix[L,Q] B0_raw[I]; 
  vector[Q] b0_int;
  vector[Q] alpha;
  vector[Q] alpha_time;
   vector[Q + choose(Q, 2)] gamma_flat; //coefficients for the (co)variances, size of S lower triangle
   vector[Q + choose(Q, 2)] gamma_flat_time; //coefficients for the (co)variances, size of S lower triangle
  real a0;
  real a0_time;
  real ran_eff[I];
  real<lower=0> ran_eff_tau;
  real<lower=0> outcome_sigma;
}

transformed parameters {
    matrix[L,Q] B0[I]; 
    matrix[L,L] Sigma[Q];
    matrix[Q, L_S] Theta[I];
    matrix[Q,Q] S[N]; 
   // matrix[L_S, L_S] Lmat[Q];
   matrix[Q,Q] Sigma0;
   matrix[L_S,L_S] theta_Sigma[Q];
   //real ran_eff[I];
    matrix[Q,Q] gamma;
   matrix[Q,Q] gamma_time;
    
   gamma = to_triangular(gamma_flat, Q);
   gamma_time = to_triangular(gamma_flat_time, Q);
  Sigma0 =  diag_pre_multiply(tau, Omega) * diag_pre_multiply(tau, Omega)';
   for(q in 1:Q){
    Sigma[q]=  diag_pre_multiply(S_tau[q], S_Omega[q]) * diag_pre_multiply(S_tau[q], S_Omega[q])';
    theta_Sigma[q]=diag_pre_multiply(theta_tau[q], theta_Omega[q]) * diag_pre_multiply(theta_tau[q], theta_Omega[q])';
   }
   for(i in 1:I){
  //   ran_eff[i] = ran_eff_tau*ran_eff_raw[i];
     for(q in 1:Q){
       B0[i][,q] = beta[q] + (diag_pre_multiply(S_tau[q], S_Omega[q])* B0_raw[i][,q]);
       Theta[i][q,] =to_row_vector(S_beta[q] + diag_pre_multiply(theta_tau[q], theta_Omega[q])*Theta_raw[i][,q]);
       }
     }
    for(n in 1:N){
       // print(Theta[ids[n]]*to_vector(S_b_spline[n,]));
        S[n] =(exp(Theta[ids[n]]*S_b_spline[n,]'))*(exp(Theta[ids[n]]*S_b_spline[n,]'))' + Sigma0;
     }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
   //outcome submodel
  a0 ~ normal(0,10);
  a0_time ~ normal(0,10);
  gamma_flat ~ normal(0,10);
  gamma_flat_time ~ normal(0,10);
 alpha ~ normal(0,10);
  alpha_time ~ normal(0,10);
  outcome_sigma ~ cauchy(0, 2.5);
  ran_eff_tau ~ cauchy(0,2.5); 

   b0_int ~ normal(0,10);
   Omega ~ lkj_corr_cholesky(1);
   tau ~ cauchy(0,1);
    for(q in 1:Q){ 
      S_beta[q] ~ normal(0, 10);
      beta[q] ~ normal(0,10);
      S_tau[q] ~ cauchy(0,1);
      S_Omega[q] ~ lkj_corr_cholesky(1);
      theta_tau[q] ~ cauchy(0, 1);
      theta_Omega[q] ~ lkj_corr_cholesky(1);
     }

    for(i in 1:I) {
      ran_eff[i] ~ normal(0, ran_eff_tau);
      for(q in 1:Q){
        B0_raw[i][,q] ~std_normal();
        Theta_raw[i][,q] ~std_normal();
      }
    }

   if(simulate == 0){ 
   vector[Q] x_mu[N]; // for each person, vector of two hormone means
   vector[N] y_mu; // for each person, vector of two hormone means
    for(n in 1:N){
      for(q in 1:Q){
        x_mu[n,q] = b0_int[q] + dot_product(B0[ids[n]][,q], b_spline[n,]);
      }
      x[n] ~ multi_normal(x_mu[n], S[n]);
      // y_mu[n] = a0+ time[n]*a0_time  + dot_product(x_mu[n,],alpha) + ran_eff[ids[n]];
      y_mu[n] = a0 + dot_product(x_mu[n,],alpha)+ time[n]*a0_time + sum(S[n] .* gamma) + sum(S[n] .* gamma_time)*time[n]  + dot_product(x_mu[n,],alpha_time)*time[n]  + ran_eff[ids[n]];
      y[n] ~ normal(y_mu[n], outcome_sigma);
    }
  }
}
// 
//   generated quantities {
//       matrix[N,Q] sim_x; 
//       vector[Q] x_mu[N];
//       vector[N] sim_y;
// //    //   //   simulate hormone data instead of using the data to generate it
//      for(n in 1:N){
//        for(q in 1:Q){
//           x_mu[n][q] = b0_int[q] + dot_product(B0[ids[n]][,q], b_spline[n,]);
//        }
//     sim_x[n,] = to_row_vector(multi_normal_rng(x_mu[n], S[n]));
//     sim_y[n] =normal_rng(a0 + dot_product(x_mu[n,],alpha)+ time[n]*a0_time + sum(S[n] .* gamma) + sum(S[n] .* gamma_time)*time[n]  + dot_product(x_mu[n,],alpha_time)*time[n]  + ran_eff[ids[n]], outcome_sigma);
//    }
// }
