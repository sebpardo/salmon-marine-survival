### Running Stan models, both hierarchical and not (simple)
### Both of which have non-centered prameterizations to improve
### sampling independence

library(tidyverse)
library(rstan)
library(bayesplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Loading data
alldatastan <- readRDS("data/alldatastan.rds")

# Loading and running hierarchical model with non-centered parameterization
modhierncp <- stan_model("analysis/survival_hier_ncp.stan")
fithierncp <- sampling(object = modhierncp, data = alldatastan, iter = 500, chains = 3,
                       control = list(adapt_delta = 0.95, max_treedepth = 15))
saveRDS(fithierncp, "data/stanfit_hier_ncp.rds")
# fithierncp <- readRDS("data/stanfit_hier_ncp.rds")


# Loading and running hierarchical model with non-centered parameterization AND new annual cv values for 1SW & 2SW separately
modhierncp_acv <- stan_model("analysis/survival_hier_ncp_annualcvs.stan")
fithierncp_acv <- sampling(object = modhierncp_acv, data = alldatastan, iter = 5000, chains = 3,
                       control = list(adapt_delta = 0.95, max_treedepth = 15))
print(fithierncp_acv)
print(fithierncp)

# fithierncpold <- fithierncp
# fithierncp <- fithierncp_acv

saveRDS(fithierncp_acv, "data/stanfit_hier_ncp_acv_full.rds")
saveRDS(fithierncpold, "data/stanfit_hier_ncp_500.rds")

saveRDS(fithierncp, "data/stanfit_hier_ncp.rds")
# fithierncp <- readRDS("data/stanfit_hier_ncp.rds")




# Investigating model object
print(fithierncp, pars = "S1")
print(fithierncp, pars = paste0("S1[",101:186,"]"))

mcmc_intervals(fithierncp_acv, pars = paste0("Pr_mu[",1:7,"]")) + 
  scale_y_discrete(labels = alldatastan$river_name)
traceplot(fithierncp_acv, pars = paste0("S1[",170:186,"]"), inc_warmup = TRUE)


# Fitting non-hierarchical model (still with ncp)
modsimplencp <- stan_model("analysis/survival_simple_ncp.stan")
fitsimplencp <- sampling(object = modsimplencp, data = alldatastan, iter = 1000, chains = 3, 
                         control = list(adapt_delta = 0.95, max_treedepth = 15))
saveRDS(fitsimplencp, "data/stanfit_simple_ncp.rds")
#fitopt <- readRDS("data/stanfit_simple_ncp.rds")

############## RUN FILE SHOULD END HERE
