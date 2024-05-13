library(rstan)
library(MASS)
library(ggplot2)
library(parallel)
options(mc.cores = parallel::detectCores(logical= FALSE))


seed = 5132024
set.seed(seed)


I = 10 ## number of observed trajectories 
L = 2 ## number of latent trajectories
J = 20 # number of observations per individual 


model_dir <- "U:/Documents/repos/blca/"

timepoints <- seq(1, J, by=1)
time <- rep(timepoints, I)
ids <- rep(seq(1, I), each=J)

id_time_df <- data.frame(cbind(ids, time))

id_time = data.frame(id_time)
ids = id_time[,2]
timepoints = id_time[,1]

### set the loadings matrix (Theta) - this is I x K, where K = size of the basis function 
## loadings matrix: 
theta1 <- c(1, 0.00, 0.25, 0.00, 0.80, 0.00, 0.50, 0.13, 0.00, 0.00)
theta2 <- c(0.00,1, 0.25, 0.40, 0.00, 0.50, 0.00, 0.01, 0.5, -0.30)
theta2 <- c(0.00,1, 0.25, 0.40, 0.00, 0.50, 0.00, 0.01, 0.5, -0.30)
theta2 <- c(0.00,1, 0.25, 0.40, 0.00, 0.50, 0.00, 0.01, 0.5, -0.30)

Theta <-cbind(theta1, theta2) # the loading matrix

### now set up the basis function for the means and the covariances 

## linear trend: 
linear_basis <- splines::bs(time[1:J], df=1, knots=c(5, 7.5, 12.5, 15))
cov_basis <- splines::bs(time[1:J], df=1, knots=c(10)) 



## coefs for trend 1: 
b1 <- rnorm(ncol(linear_basis), 0, 1)
b2 <- rnorm(ncol(linear_basis), 0, 0.25)


B = cbind(b1, b2)

l1 <- data.frame(linear_basis)
### multiply the basis function by Theta now:

Lambda = Theta%*%t(B)
lin
x_mu = t(sapply(seq_along(timepoints), function(i){Lambda%*%t(l1[i,])}))


### plot these trajectories: 
latent_df <- data.frame(cbind(time=seq(1,J)), x_mu)
latent_melt <- reshape2::melt(latent_df, id.vars=c("time"), variable.name="id", value.name="trajectory")

## this is plotting the means of the trajectories 
ggplot(data=latent_melt, aes(x=time, y=trajectory,group=id, color=id)) + 
  geom_line(linewidth=1.5) + 
  geom_point(size=3) + 
  facet_wrap(~id)

## generate now the covariance regression 
## we'll just do a linear trend with 1 knot 
# cov_basis <- splines::bs(time[1:J], df=1, knots=c(10))
# ## coefs for trend 1: 
# b1_cov <- rnorm(ncol(cov_basis), 0, 0.5)
# b2_cov  <- rnorm(ncol(cov_basis), 0, 0.05)

### create Sigma0 
sigma0 <- rlnorm(I, 0, 0.05)
Sigma0 <- diag(sigma0)

Sigma = Lambda%*%t(Lambda) + Sigma0


## simulate the means of the markers based on the loadings matrix and the latent traits  
sim_mu <- t(sapply(seq(nrow(scale_time_df)), function(i){
  return(Lambda%*%eta_i[i,])
}))

sigma_q <- exp(rnorm(Q,0, 1))
Sigma_q <- diag(sigma_q)


sim_markers <- t(sapply(seq_along(ids), function(i){MASS::mvrnorm(1, sim_mu[i,], Sigma_q)}))

## compile model 
compiled_model <- stan_model(paste0(model_dir, "factor_model.stan"))

## set results directory
results_dir <- "U:/Documents/repos/blca/"

## sample from model 
sim_out <- sampling(compiled_model,
                    # include = TRUE,
                    sample_file=paste0(results_dir, 'rep0_',seed, "_", Sys.Date(), '_model_samples.csv'), #writes the samples to CSV file
                    iter =2000,
                    warmup=1000, #BURN IN
                    chains = 2,
                    seed = seed,
                    control = list(max_treedepth = 30,
                                   adapt_delta=0.95),
                    data = list(Q = Q,
                                X = sim_markers,
                                N = N,
                                K = K, 
                                I = N,
                                P = P,
                                Z= scale_time_df,
                                # Psi = Psi,
                                Sigma = Sigma_q,
                               # eta = t(eta_i),
                                ids =ids,
                               simulate = 0
                    ))

# 
# 
# sample_file_names <-c("rep0_30_2024-04-17_model_samples_1", "rep0_30_2024-04-17_model_samples_2")
# sim_out <- read_stan_csv(paste0(results_dir, sample_file_names,".csv"))
# 
# rstan::traceplot(sim_out, pars=c("lambda"))
# sim_summary <- data.frame(summary(sim_out, pars=c("lambda", "Lambda"))$summary)
