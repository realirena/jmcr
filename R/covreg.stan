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
  int<lower=1> L; // number of latent factors 
  int<lower=1> K; // dim of the basis expansion on time 
  int<lower=1> id[N]; //array of subject ids (length N in order to do the longitudinal estimation)
  matrix[N,I] x; // the variable to be forecasted - in this case, number of applicants 
  vector[N] time; // vector that stores the time index 
  matrix[N, K] mu_design; // this is the design matrix for how the means of the factors evolve 
  matrix[N, K] S_design; // this is the design matrix for how the factor variances and correlations evolve
}

// The parameters accepted by the model.
parameters {
  vector[L] S_beta; 
  vector[L] Theta_raw;  // prior correlation (for the covariance matrix of Bijs)
  vector<lower=0>[L] S_tau;
  cholesky_factor_corr[I] Omega;  // prior correlation (for the covariance matrix of Bijs)
  vector<lower=0>[I] tau;
  cholesky_factor_corr[L] S_Omega;  // prior correlation (for the covariance matrix of Bijs)
  vector<lower=0>[L] theta_tau; 
  cholesky_factor_corr[L] theta_Omega;  // prior correlation (for the covariance matrix of Bijs)
}

transformed parameters {
   vector[L] B0; 
   matrix[L,L] Sigma;
   matrix[I, L] Theta;
   matrix[I,I] S[N]; 
   matrix[I,I] Sigma0;
   matrix[L,L] theta_Sigma;
   matrix[I, K] Lambda[N];
   
   Sigma0 =  diag_pre_multiply(tau, Omega) * diag_pre_multiply(tau, Omega)';
   theta_Sigma=diag_pre_multiply(theta_tau, theta_Omega) * diag_pre_multiply(theta_tau, theta_Omega)';
   Sigma=  diag_pre_multiply(S_tau, S_Omega) * diag_pre_multiply(S_tau, S_Omega)'; // compute the covariance matrix for the regression coefs 
   
  for(i in 1:I){
    Theta[i,] =to_row_vector(S_beta + diag_pre_multiply(theta_tau, theta_Omega)*Theta_raw); // this gets the factor loadings matrix 
  
  } 
   
  for(n in 1:N){
    Lambda[n] =exp(Theta*S_design[n,]');
    S[n] =Lambda[n]*Lambda[n]' + Sigma0; // computes the covariance matrix at the n'th timepoint 
  }
}

model {
    Omega ~ lkj_corr_cholesky(1);
    tau ~ cauchy(0,1);
    S_beta ~ normal(0, 10);
    S_tau ~ cauchy(0,1);
    S_Omega ~ lkj_corr_cholesky(1);
    theta_tau ~ cauchy(0, 1);
    theta_Omega ~ lkj_corr_cholesky(1);
    Theta_raw ~std_normal();
    for(n in 1:N){
      x[n] ~ multi_normal((Lambda[n]*mu_design[n,]')', S[n]);
    }
  }


// 
//   generated quantities {
//       matrix[N,I] sim_x; 
//       vector[I] x_mu[N];
//       vector[N] sim_y;
// //    //   //   simulate hormone data instead of using the data to generate it
//      for(n in 1:N){
//        for(q in 1:I){
//           x_mu[n][q] = b0_int[q] + dot_product(B0[ids[n]][,q], mu_design[n,]);
//        }
//     sim_x[n,] = to_row_vector(multi_normal_rng(x_mu[n], S[n]));
//     sim_y[n] =normal_rng(a0 + dot_product(x_mu[n,],alpha)+ time[n]*a0_time + sum(S[n] .* gamma) + sum(S[n] .* gamma_time)*time[n]  + dot_product(x_mu[n,],alpha_time)*time[n]  + ran_eff[ids[n]], outcome_sigma);
//    }
// }
