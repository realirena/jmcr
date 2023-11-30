rm(list=ls())
library(rstan)
library(dplyr)
library(reshape2)
library(splines)
library(MASS)
library(LaplacesDemon)
library(Matrix)
library(nimble)
#library(matrixcalc)
options(mc.cores = parallel::detectCores(logical= FALSE))
# rstan_options(auto_write = TRUE)
taskid <- 2022
## set up the data simulation parameters:
seed = 3*taskid + 10

Q = 2 # no of predictor markers 
I = 100 # no of subjects
J = 5 # no of time points (days)

b0_int_est <- data.frame(read.csv("/b0_int_est_0627.csv"))
## samples for the means and v-cov of the B-spline 
theta_Sigma_est <- data.frame(read.csv("/beta_Sigma_est_0627.csv"))
Sigma_est <- data.frame(read.csv("/Sigma_est_0627.csv"))
beta_est <- data.frame(read.csv("/beta_est_0627.csv"))
Sigma0_est <- data.frame(read.csv("/Sigma0_est_0627.csv"))
beta1 <-  beta_est[1:5, 2]
beta2 <-  beta_est[6:10, 2]

beta <- list()
beta[[1]] <- beta1
beta[[2]] <-  beta2

Sigma1 <-  matrix(Sigma_est[1:25,2], ncol=5, nrow=5)
Sigma2 <-  matrix(Sigma_est[26:50,2], ncol=5, nrow=5)
Sigma <- list()
Sigma[[1]] <- Sigma1
Sigma[[2]] <- Sigma2

b0_int <- b0_int_est[1:2,2]


Ti = rep(J, I)
#Ti = sample(1:J, I, replace=T) ##no visits per subject
## placeholder for the simulated time to FMPs
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
for(i in 1:I){
  time_fmp_i = seq(-J,J, by=1)
  b = min(as.numeric(rpois(1, lambda = 5)),2*J-1)
  if(i==1){
    if(b==0){
      id_time = time_fmp_i
    }else{
      id_time = time_fmp_i[-sample(1:length(time_fmp_i), b)]
    }
    scale_time = range01(id_time)
    scale_time_df <- cbind(1, scale_time)
    head(scale_time_df)
    id_time = cbind(id_time, rep(i, length(id_time)))
  } else{
    if(b==0){
      id_time_tmp = time_fmp_i
    }else{
      id_time_tmp = time_fmp_i[-sample(1:length(time_fmp_i), b)]
    }
    scale_time_tmp <- cbind(1, range01(id_time_tmp))
    id_time_tmp <- cbind(id_time_tmp, rep(i, length(id_time_tmp)))
    id_time = rbind(id_time, id_time_tmp)
    scale_time_df = rbind(scale_time_df, scale_time_tmp)
  }
}
Ti = tabulate(id_time[,2])

N= nrow(id_time)

id_time = data.frame(id_time)
ids = id_time[,2]
timepoints = id_time[,1]

## place knots at -7, -2, 2, and 7 years to FMP 
knts = c(-7, -2,2, 7)
num_knots = length(knts)
degree= 1
num_basis = num_knots + degree 
## generate the individual mean coefficients 
B = list()
for(i in 1:I){
  mat <- matrix(nrow=num_basis, ncol=Q)
  for(q in 1:Q){
    mat[,q] <- MASS::mvrnorm(1,beta[[q]], Sigma[[q]])
    # mat[,q] = abs(mat[q])
  }
  B[[i]] <- mat
}


## simulate the mean trajectories based on the timepoints
for(i in 1:I){
  if(i==1){
    time_fmp =id_time[id_time$V2==i,]$id_time
    b_spline  <-splines::bs(time_fmp, degree=degree, knots=knts) ## remove if we add an intercept term 
    sim_mu = matrix(nrow=Ti[i], ncol=Q)
    for(q in 1:Q){
      sim_mu[,q] = b0_int[q] +b_spline%*%B[[i]][,q]
    }
    f_time_fmp = cbind(1, time_fmp)
  } else {
    time_fmp_i =id_time[id_time$V2==i,]$id_time
    b_spline_tmp <-splines::bs(time_fmp_i, degree=degree, knots=knts)## remove if we add an intercept term
    sim_mu_tmp = matrix(nrow=Ti[i], ncol=Q)
    for(q in 1:Q){
      sim_mu_tmp[,q] = b0_int[q] +b_spline_tmp%*%B[[i]][,q]
    }
    b_spline = rbind(b_spline, b_spline_tmp)
    sim_mu = rbind(sim_mu, sim_mu_tmp)
    f_time_fmp = rbind(f_time_fmp, cbind(1, time_fmp_i))
  }
}

S_beta = list()
S_beta[[1]] = c(-1.5, -5)
S_beta[[2]] = c(-1,-0.85)
theta_Sigma <- list()

## preset true values of the cov mat for the individual covariances 
theta_Sigma[[1]]<- matrix(c(1, 0.002, 0.002, 0.001), ncol=2)
theta_Sigma[[2]]<-  matrix(c(1.24, 0.002, 0.002, 0.07), ncol=2)


Sigma0 = matrix(c(0.3, -0.2, -0.2, 0.6), ncol=Q)

f_time_fmp <- data.frame(cbind(1, id_time$id_time))
f_time_fmp$id <- ids

## generate the variance and covariance time-trajectories 
S_knts = 0
S_degree= 1
P_S = length(S_knts) + S_degree + 1 
for(i in 1:I){
  if(i==1){
    time_fmp =f_time_fmp[f_time_fmp$id==i,2]
    S_b_spline  <-splines::bs(time_fmp, degree=S_degree, knots=S_knts) ## remove if we add an intercept term
  } else {
    time_fmp_i = f_time_fmp[f_time_fmp$id==i,2]
    S_b_spline_tmp <-splines::bs(time_fmp_i, degree=S_degree, knots=S_knts)## remove if we add an intercept term
    S_b_spline = rbind(S_b_spline, S_b_spline_tmp)
  }
}

## function to create the individual time-varying covariance matrices 
S_b_spline <- cbind(1,f_time_fmp[,2])
P_S = ncol(S_b_spline)
get_tv_vcov = function(S_b_spline, Q,  P_S, ids, N, theta_Sigma, Sigma0){
  S = list()
  Theta <- list()
  for(i in 1:I){
    theta <- matrix(NA, nrow=Q, ncol=P_S)
    for(q in 1:Q){
      theta[q,] = MASS::mvrnorm(1, S_beta[[q]],theta_Sigma[[q]])
    #  theta[q,] = abs(theta[q,])
    }
    Theta[[i]] <- theta
  }
  
  for(n in 1:N){
    xi <- S_b_spline[n,]
    S[[n]]=exp(Theta[[ids[n]]]%*%xi)%*%t(exp(Theta[[ids[n]]]%*%xi)) + Sigma0
  }
  return(S)
}

## apply cov function to each individual 
S <-  get_tv_vcov(S_b_spline = S_b_spline,
                  Q = Q, P_S = P_S,
                  Sigma0 = Sigma0,
                  theta_Sigma = theta_Sigma,
                  N=length(ids), ids)


### generate the hormone data (time-varying means and variances)
sim_x = sapply(seq_along(ids), function(i){MASS::mvrnorm(1, sim_mu[i,],S[[i]])})
S_lower <- sapply(seq_along(ids), function(i){
  mat = S[[i]]
  mat[upper.tri(mat)==T] <- 0
  return(mat)
})
## unroll S into a design matrix: 
S_design <- matrix(unlist(S_lower), ncol=Q*Q, byrow=T)
S_design <- S_design[,-3]

a0 = -0.5
alpha <- c(-1.25, 1.5)
a0_time <- c(0.9)
gamma <- c(-0.75, 1.5, -0.25)
alpha_time <- c(0.1, -0.03)
gamma_time <- c(-0.05, 0.1,  -0.08)
## add a random effect for the outcome
ran_eff_tau = 0.25
ran_eff = sapply(seq_along(1:I), function(i){rnorm(1, 0, ran_eff_tau)})

## set the variance of the outcome:
outcome_sigma <- 0.10
### generate the outcome dat

sim_mu_outcome <- sapply(seq_along(ids), function(i){a0 + a0_time*f_time_fmp[i,2]  + sim_mu[i,]%*%alpha + S_design[i,] %*%gamma+ (S_design[i,]%*%gamma_time)*f_time_fmp[i,2] + (sim_mu[i,]%*%alpha_time)*f_time_fmp[i,2]+ ran_eff[ids[i]]})
#sim_mu_outcome <- sapply(seq_along(ids), function(i){a0 + a0_time*time[i]  + sim_mu[i,]%*%alpha+ ran_eff[ids[i]]})

sim_outcome <- sapply(seq_along(ids), function(i){rnorm(1, sim_mu_outcome[i],outcome_sigma)})

## compile model 
compiled_model <- stan_model("/tvv_joint.stan")

## pass in initial values (unsure if this actually speeds up anything)
initf1 <- function() {
  list(beta =beta,
       Sigma = Sigma,
       Sigma0 = Sigma0,
       theta_Sigma = theta_Sigma,
       B0 = B,
       S_beta = S_beta,
       S = S,
       a0=a0,
     alpha = alpha,
      alpha_time = alpha_time,
      a0_time = a0_time,
      ran_eff_tau = ran_eff_tau,
      gamma_flat = gamma,
      gamma_flat_time = gamma_time,
       outcome_sigma = outcome_sigma
  )
}

## set results directory
results_dir <- "/insertpathhere/"

## sample from model 
sim_out <- sampling(compiled_model,
                     # include = TRUE,
                     sample_file=paste0(results_dir, 'rep_',taskid, '_model_samples.csv'), #writes the samples to CSV file
                     iter =2000,
                     warmup=1000, #BURN IN
                     chains = 2,
                     seed = seed,
		    init = initf1,
                     control = list(max_treedepth = 30,
                                    adapt_delta=0.95),
                     data = list(Q = Q,
                                 I = I,
                                 N = N,
                                 L= num_basis,
                                 L_S = ncol(S_b_spline),
                                 ids=ids,
                                  Sigma0 = Sigma0,
                                 b_spline=b_spline,
                                 S_b_spline= S_b_spline,
                                 x = t(sim_x),
				                         time = f_time_fmp[,2],
				                         y = sim_outcome,
                                 simulate = 0
                     ))


