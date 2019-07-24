# Convert Scrogster's JAGS model to greta model

# NG notes: for loops and indexing are slower in greta than vectorising things
# and using matrix algebra, because it prevents tensorflow form parallelising
# operations. They also increase the overhead in running the model. So a lot of
# things here are pared down to make it more efficient in processing.

# NG cont: The for loop for the AR1 process in particular looks very different,
# beecause everything that could be pre-computed has been.

# NG: I also replaced the indicator variables on relevant lags with a weighted
# sum, which should mix much better.

# NG: I dropped the foxes, and the observation model is Poisson (lognormal), not
# zero-inflated Poisson (lognormal), though we could potentially do the ZIP
# model

# alternative way to simulate AR1. Can expand the iterative equation:
#   X_t = \rho X_{t-1} + epsilon_i
# to
#   X_t = \Sum_{i=0}^{n}(\epsilon_{t - n} \rho ^ n)
# which can be simulated with elementwise matrix operations, and a matrix
# multiplication. Sigma can be a vector with a different element per site
ar1 <- function (rho, n_times, n_sites = 1, sigma = 1,
                 innovations = normal(0, 1, dim = c(n_times, n_sites))) {
  
  # matrix of time modifiers
  t_seq <- seq_len(n_times)
  t_mat <- outer(t_seq, t_seq, FUN = "-")
  t_mat <- pmax(t_mat, 0)
  
  # which elements to include (don't want upper triangular ones)
  mask <- lower.tri(t_mat, diag = TRUE)
  
  # matrix of rho ^ n contributions
  rho_mat <- (rho ^ t_mat) * mask
  
  # multiply by scaled innovations to get timeseries
  if (length(sigma) == n_sites) {
    innovations_scaled <- sweep(innovations, 2, sigma, "*")
  } else if(length(sigma) == 1) {
    innovations_scaled <- innovations * sigma
  } else {
    stop ("sigma must have either length 1 or length n_sites")
  }
  rho_mat %*% innovations_scaled
  
}

load("prepped_data.Rdata")

# make / adjust data object
n_obs <- hier_dat$Nobs
site_code <- hier_dat$site.code
n_sites <- length(unique(hier_dat$site.code))
n_lag <- dim(rain_lag_array)[3]
n_times <- 40

library(greta)

# helper functions for Scroggie's positive and continuous priors.

# Scroggie specified these as smoothed spike and slab distributions (a Student
# distribution with 1 df, so that they have fat tails to allocate more weight of
# probability to extreme values).  We couldn't get good mixing with those, so
# have retreated to Normal distributions.  Could go for something intermediate I
# guess (a Student distribution with higher df), but as NG points out you have
# to stop and think about whether in principle you believe those parameters
# could take extreme values.

continuous <- function(sd = 2.5, dim = NULL) {
  normal(0, sd, dim = dim)
}

positive <- function(sd = 1, dim = NULL) {
  normal(0, sd, truncation = c(0, Inf), dim = dim)
}

# random site effects on r.max for rabbit
r_mean_rabbits <- continuous()
site_sd_rabbits <- positive()

# could use this:
#   site_r_effect_rabbits <- normal(r_mean_rabbits, site_sd_rabbits, dim = n_sites)
# but instead use decentred version of this hierarchical structure, for more efficient samplers
# keep centred versions for plotting, or whatever
site_r_effect_rabbits_raw <- normal(0, 1, dim = n_sites)
site_r_effect_rabbits_centred <- site_r_effect_rabbits_raw * site_sd_rabbits
site_r_effect_rabbits <- r_mean_rabbits + site_r_effect_rabbits_centred

# lagged rain effect ----
# create a rain effect matrix by sweeping a vector of probabilities over the
# rain array

# create a simplex of weights on different rain lags
lag_weights_raw <- uniform(0, 1, dim = n_lag)
lag_weights <- lag_weights_raw / sum(lag_weights_raw)

# get a weighted rain effect matrix, by applying these weights to the different
# lags, and summing them. This is the same as marginalising the discrete lags
# analytically, or assuming that all lags have some effect, and weighting them
# all probabilistically. Should mix better than the discrete version and be as
# interpretable
rain_lag_array <- rain_lag_array[, seq_len(n_times), ]
rain_lag_array[is.na(rain_lag_array)] <- 0
# reshape the array to a matrix, do a matrix multiply with thge weights, and then reshape :O
rain_lag_array_long <- rain_lag_array
dim(rain_lag_array_long) <- c(n_sites * n_times, n_lag)
weighted_rain_lags_long <- rain_lag_array_long %*% lag_weights
weighted_rain_lags <- weighted_rain_lags_long
dim(weighted_rain_lags) <- c(n_sites, n_times)

# phew!

# combine with the rain coefficient to get the rain effect ----
rain_coef <- continuous()
rain_effect <- weighted_rain_lags / 10 * rain_coef


# get the temporal effects ---- 

# looks like there was a bug in defining these variables in prepped_data!
# guessing that they should look like this instead:
postrip <- c(rep(0, 13), rep(1, 27))
winter <- rep(c(0, 1), 20)

winter_coef <- continuous()
postrip_coef <- continuous()

winter_effect <- winter * winter_coef
postrip_effect <- postrip * winter_coef

temporal_effect <- winter_effect + postrip_effect

# build non-dynamic part of rabbit change matrix
env_matrix_rabbits <- sweep(rain_effect, 2, temporal_effect, "+")
non_dynamic_rabbits_t <- sweep(env_matrix_rabbits, 1, site_r_effect_rabbits, "+")

# use an AR1 model for the rabbit dynamics
auto_coef <- continuous()
proc_sd_rabbits <- positive()

# create a matrix of standard normal temporal deviates here, for decentring and
# efficiency
temporal_deviates_raw <- normal(0, 1, dim = c(n_times - 1, n_sites))
temporal_deviates  <- temporal_deviates_raw * proc_sd_rabbits

# combine these with the other covariates of the (density-independent part of
# the) log growth rates
log_growth_rates <- t(non_dynamic_rabbits_t)[-1, ] + temporal_deviates

# initial population sizes
log_mu_rabbits_init_raw <- normal(0, 1, dim = c(1, n_sites))
log_mu_rabbits_init <- log(4) + log_mu_rabbits_init_raw * 0.25

# solve the AR(1) process using these density-independent bits and the density
# dependent part (autoregressive component controlled by auto_coef)
innovations <- rbind(log_mu_rabbits_init,
                     log_growth_rates)

log_mu_rabbits <- ar1(
  rho = auto_coef,
  n_times = n_times,
  n_sites = n_sites,
  sigma = 1,
  innovations = innovations
)

# pull the relevant bits out of the matrix of site/time combinations, into long
# form matching the observations
indices <- cbind(hier_dat$obs_time, hier_dat$site.code)
log_mu_rabbits_obs <- log_mu_rabbits[indices]

# observation models
log_lambda_rabbits <- log_mu_rabbits_obs +  log(hier_dat$trans.length / 1000)#  + surv_err_rabbit
expected_rabbits <- exp(log_lambda_rabbits)

# poisson
# distribution(hier_dat$rabbit.count) <- poisson(expected_rabbits)

# convert to negative binomial observation parameters
size <- normal(0, 10, truncation = c(0, Inf))
prob <- 1 / (1 + expected_rabbits / size)
distribution(hier_dat$rabbit.count) <- negative_binomial(size, prob)

set.seed(2019-07-03)

m <- model(winter_coef, rain_coef, postrip_coef, auto_coef)
system.time(
  draws <- mcmc(m, sampler = hmc(Lmin = 30, Lmax = 40), warmup = 1000, chains = 32)
)
plot(draws)
