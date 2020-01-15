library(tidyverse)
library(rstan)
library(bayesplot)

# loading data
alldatastan <- readRDS("data/alldatastan.rds")

# loading stanfit objects
fithierncp <- readRDS("data/stanfit_hier_ncp.rds") # hierarchical
#fitsimplencp <- readRDS("data/stanfit_simple_ncp.rds")   # non-hierarchical

summary(fithierncp)$summary[,"Rhat"]
range(summary(fithierncp)$summary[,"Rhat"])

neffs <- summary(fithierncp)$summary[,"n_eff"]
neffs[neffs < 1000]

range(summary(fithierncp)$summary[,"n_eff"])

print(fithierncp, pars = "Pr_mu")
print(fithierncp, pars = "logitPr_sigma")

print(fithierncp, pars = "Z1")

ratios_ncp <- neff_ratio(fithierncp)
mcmc_neff(ratios_ncp, size = 2)

500/7500
7