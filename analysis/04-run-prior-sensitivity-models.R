# Running Stan models with differing prior assumptions to assess sensitivity
# to prior choice
#  Sensitivity tests: 
# (1) use a relatively weak informative prior for Z (equations 12 and 13), e.g., using sigma = 0.5 or 0.8; 
# (2) increase/decrease the median values for Z by 10% or 20%; 
# (3) you may also do the similar tests for parameter Pg,r,t in equations 14-16. 
  
library(tidyverse)
library(rstan)
library(bayesplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Loading data
alldatastan <- readRDS("data/alldatastan.rds")

# (1) use a relatively weak informative prior for Z (equations 12 and 13), e.g., using sigma = 0.5 
# (2) increase/decrease the median values for Z by 10% or 20%; 
# lower values of mu (towards -Inf result in higher S values skewed towards 1)
# values closer to zero center S around 0.5. 
# Basically, this is weakening the prior towards a median of 0.5 with wider tails
alldatastan_zmu_h <- alldatastan
alldatastan_zmu_h$z1_mu <- 0.5
alldatastan_zmu_h$z1_sd <- 0.5
alldatastan_zmu_h$z2_mu <- 0.1
alldatastan_zmu_h$z2_sd <- 0.5

x <- seq (0.001, 0.999, by =0.001)
hist(exp(-rlnorm(5000, alldatastan_zmu_h$z1_mu, alldatastan_zmu_h$z1_sd)))
hist(exp(-rlnorm(5000, alldatastan_zmu_h$z2_mu, alldatastan_zmu_h$z2_sd)))

# In the original stan model:
# Z1 ~ lognormal(1, 0.22);    
# Z2 ~ lognormal(0.2, 0.3);

plot(x, dlnorm(-log(x), 1, 0.22))
plot(x, dlnorm(-log(x), alldatastan_zmu_h$z1_mu, alldatastan_zmu_h$z1_sd))
plot(x, dlnorm(-log(x), 0.2, 0.3))
plot(x, dlnorm(-log(x), alldatastan_zmu_h$z2_mu, alldatastan_zmu_h$z2_sd))


modhierncp_acv_priors <- stan_model("analysis/survival_hier_ncp_annualcvs_datapriors.stan")

fithierncp_zsdmu <- sampling(object = modhierncp_acv_priors, 
                           data = alldatastan_zsd, iter = 5000, chains = 3,
                           control = list(adapt_delta = 0.95, max_treedepth = 15))
print(fithierncp_zsdmu)

saveRDS(fithierncp_zsdmu, "data/stanfit_hier_ncp_acv_zsdmu.rds")
fithierncp_zsdmu <- fithierncp_zsd



# (3) Change population-level priors of mu_Pg,r and sigma_Pg,r to all be like
# non-1SW dominated prior probit(0, 2.8)

alldatastan_pr_weak <- alldatastan

(alldatastan_pr_weak$logis_mu <- alldatastan$logis_mu * 0.6956522)
(alldatastan_pr_weak$logis_sigma <- alldatastan$logis_sigma + 0.4)

alldatastan_pr_weak$logis_mu <- rep(alldatastan$logis_mu[1], 7)
alldatastan_pr_weak$logis_sigma <-  rep(alldatastan$logis_sigma[1], 7)

#plot(x, dlogis(qnorm(x), 2.3, 0.4))
#plot(x, dlogis(qnorm(x), 1.6, 0.8))

plot(x, dlogis(qnorm(x), 0, 2.8)) 
#plot(x, dlogis(qnorm(x), 0, 3.2)) 

# plot(x, dlogis(qlogis(x), 2.3, 0.4))
# plot(x, dlogis(qnorm(x), 1.7, 0.4))
# plot(x, dlogis(qnorm(x), 2.3, 0.4))
# 
# plot(x, dlogis(qlogis(x), 0, 0.8)) 
# plot(x, dlogis(qnorm(x), 0, 0.5)) 
# plot(x, dlogis(qnorm(x), 0, 0.8)) 


# Z1 and Z2 priors used in main model
alldatastan_pr_weak$z1_mu <- 1
alldatastan_pr_weak$z1_sd <- 0.22
alldatastan_pr_weak$z2_mu <- 0.2
alldatastan_pr_weak$z2_sd <- 0.3


fithierncp_prweak <- sampling(object = modhierncp_acv_priors, 
                             data = alldatastan_pr_weak, iter = 5000, chains = 3,
                             control = list(adapt_delta = 0.95, max_treedepth = 15))


print(fithierncp_prweak)
saveRDS(fithierncp_prweak, "data/stanfit_hier_ncp_acv_prweak.rds")


