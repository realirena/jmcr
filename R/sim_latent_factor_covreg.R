rm(list=ls())
library(rstan)
library(MASS)
library(ggplot2)
library(parallel)
options(mc.cores = parallel::detectCores(logical= FALSE))

seed = 5132024
set.seed(seed)


I = 10 ## number of observed trajectories 
L = 2 ## number of latent trajectories
N = 20 # number of observations per individual 


model_dir <- "U:/Documents/repos/jmcr/R/"


### simulating timepoints 
timepoints <- seq(1, N, by=1)
time <- rep(timepoints, I)
ids <- rep(seq(1, I), each=N)


### set the loadings matrix (Theta) - this is I x K, where K = size of the basis function 

## loadings matrix: 
theta1 <- c(1, 0.00, 0.25, 0.00, 0.80, 0.00, 0.50, 0.13, 0.00, 0.00)
theta2 <- c(0.00,1, 0.25, 0.40, 0.00, 0.50, 0.00, 0.01, 0.5, -0.30)
# 
# theta3  <- c(0.00, 0.00,1, 0.3, -0.25, 0.4, 0.1, 0.15, 0.00, 0.00)
# theta4  <- c(0.00, 0.00, 0.00,1, 0.30, 0.00, 0.00, 0.00, 0.00, 0.00)

Theta <-cbind(theta1, theta2) # the loading matrix

### now set up the basis function for the means and the covariances 


### assume that the factors have the same covariance regression for now (linear time trend)
## put 1 knot at time = 10

### K = dimension of the covariance regression basis functions
K =2 #intercept and slope 

cov_basis <- list(cb1 = cbind(1, time[1:N]),
                  cb2 = cbind(1, time[1:N]))

## coefs for trend 1: 
b1_cov <- c(0.24, -0.05) ## high starting, decreasing variance over time
b2_cov <- c(-0.2, 0.16) ## low variance, bu increasing over time 
B_cov <- cbind(b1_cov, b2_cov)

## these are the latent "covariance" trajectories 
eta = matrix(nrow=N, ncol=K)
eta[,1] <- cov_basis[[1]]%*%B_cov[,1]
eta[,2] <- cov_basis[[2]]%*%B_cov[,2]

## the "loadings matrix" for the covariance regression and the mean 
Lambda <- list()
for(n in 1:N){
  tmp <- matrix(nrow=I, ncol=L)
  for(l in 1:L){
    tmp[,l]  <- exp(Theta[,l]%*%t(eta[n,l]))
  }
  Lambda[[n]] <- tmp
}

### preset Sigma0 
sigma0 <- rlnorm(I, 0, 0.05)
Sigma0 <- diag(sigma0)

### Hold the time-varying covariances: 

S <- list()
for(n in 1:N){
  S[[n]] <- Lambda[[n]]%*%t(Lambda[[n]]) + Sigma0
}

## linear trend:
linear_basis <- list()
linear_basis[[1]] <- as.matrix(splines::bs(time[1:N], degree=1, knots=c(5, 10, 15)))
linear_basis[[2]] <- as.matrix(splines::bs(time[1:N], degree=1, knots=c(5,10,12)))

## coefficients for the basis function 
b1 <- rnorm(ncol(linear_basis[[1]]), 0, 1)
b2 <-  rnorm(ncol(linear_basis[[2]]), 0, 0.5)

B <- as.matrix(cbind(b1, b2))

phi <- matrix(nrow=N, ncol=L)

for(l in 1:L){
  phi[,l] = linear_basis[[l]]%*%B[,l]
}

### multiply the basis function by Theta now:
sim_mu <- t(sapply(seq(N), function(i){
  return(Lambda[[i]]%*%phi[i,])
}))


latent_df <- data.frame(cbind(time=seq(1,N)), sim_mu)
latent_melt <- reshape2::melt(latent_df, id.vars=c("time"), variable.name="id", value.name="mean_trajectory")

## this is plotting the means of the trajectories 
ggplot(data=latent_melt, aes(x=time, y=mean_trajectory,group=id, color=id)) + 
  geom_line(linewidth=1.5) + 
  geom_point(size=3) + 
  facet_wrap(~id)


sim_x <- sapply(seq(N), function(i){
  return(MASS::mvrnorm(1, sim_mu[i,], S[[i]]))
})

latent_df <- data.frame(cbind(time=seq(1,N)), sim_x)
latent_melt <- reshape2::melt(latent_df, id.vars=c("time"), variable.name="id", value.name="trajectory")

## this is plotting the means of the trajectories 
ggplot(data=latent_melt, aes(x=time, y=trajectory,group=id, color=id)) + 
  geom_line(linewidth=1.5) + 
  geom_point(size=3) + 
  facet_wrap(~id)

## set results directory
results_dir <- "U:/Documents/repos/blca/"
compiled_model <- stan_model(paste0(model_dir, "covreg.stan"))

## sample from model 
sim_out <- sampling(compiled_model,
                    # include = TRUE,
                    #sample_file=paste0(results_dir, 'rep0_',seed, "_", Sys.Date(), '_model_samples.csv'), #writes the samples to CSV file
                    iter =200,
                    warmup=100, #BURN IN
                    chains = 1,
                    seed = seed,
                    control = list(max_treedepth = 30,
                                   adapt_delta=0.95),
                    data = list(
                                x = t(sim_x),
                                N = N,
                                I = I,
                                L = L,
                                K= 4,
                                S_K = 2,
                                mu_design=linear_basis,
                                S_design = cov_basis,
                               simulate = 0
                    ))

# 
# 
# sample_file_names <-c("rep0_30_2024-04-17_model_samples_1", "rep0_30_2024-04-17_model_samples_2")
# sim_out <- read_stan_csv(paste0(results_dir, sample_file_names,".csv"))
# 
# rstan::traceplot(sim_out, pars=c("lambda"))
# sim_summary <- data.frame(summary(sim_out, pars=c("lambda", "Lambda"))$summary)
