# Merge data from all populations into structure to input to TMB hierachical models
# as well as create object for Stan models

library(tidyverse)
library(rlist)
#library(data.table)
library(rowr)

# NEED TO CHECK FOR -Inf VALUES!!!

saintjeandata <- readRDS(file = "data/saintjeandata.rds")
trinitedata <- readRDS(file = "data/trinitedata.rds")
connedata <- readRDS(file = "data/connedata.rds")
nashwaakdata <- readRDS(file = "data/nashwaakdata.rds")
wabdata <- readRDS(file = "data/wabdata.rds")
campdata <- readRDS(file = "data/campbelltondata.rds")
lahavedata <- readRDS(file = "data/lahavedata.rds")

# Converting lists for separate population into single list with  matrices for each data type 
# where each population is a separate row
# based on answer in 
# https://stackoverflow.com/questions/57401290/bind-vectors-across-lists-to-single-list-of-matrices

alldata <- list(saintjeandata, trinitedata, connedata, nashwaakdata, wabdata, campdata, lahavedata) %>% 
  purrr::transpose() %>%
  map(~ reduce(.x, cbind.fill, fill = NA) %>% 
        t %>% 
        `row.names<-`(NULL)) %>%
  list.remove("allreturns")

saveRDS(alldata, file = "data/alldata.rds")

# Creating list for Stan models

# listnames <- names(alldata)

alldatastan <- alldata %>%
  lapply(t) %>%
  lapply(as.vector) %>%
  lapply(na.omit) %>%
  lapply(as.vector) # To remove na.omit class and subsequent join warnings

# names(alldatastan) <- listnames

# Saving filtered data frame in new object (connereturns)
alldatastan$N <- length(alldatastan$nyears)
alldatastan$K <- sum(alldatastan$nyears)

# re-adding matrix of years
alldatastan$yearmat <- alldata$year

grilse_rivers <- c("Conne River", "Western Arm Brook", "Campbellton River")

# logistic priors for population level Pr 
alldatastan$logis_mu <- ifelse(alldatastan$river_name %in% grilse_rivers, 2.3, 0)
alldatastan$logis_sigma <- ifelse(alldatastan$river_name %in% grilse_rivers, 0.4, 2.8)

alldatastan$logis_sigma
alldatastan$logsmolts_cv

lapply(alldatastan, length)

saveRDS(alldatastan, file = "data/alldatastan.rds")

names(alldatastan)

# merging returns data
qcreturns <- bind_rows(saintjeandata$allreturns, trinitedata$allreturns) %>% 
  rename(M1 = returns_grilse, M2 = returns_2sw_est, M3 = returns_3sw_est, RS = returns_rs_est) %>%
  mutate(totalreturns = M1 + M2 + M3 + RS)

nlreturns <- bind_rows(connedata$allreturns, wabdata$allreturns, campdata$allreturns) %>%
  rename(smolts = smolt, M1 = estM1, M2 = estM2p, RS = estRS)

nbreturns <- nashwaakdata$allreturns %>%
  select(year, smolts = smolt_mode, M1 = est_1SW, M2 = est_2SW, RS = est_RS, river_name)  %>%
  mutate(totalreturns = M1 + M2 + RS)

nsreturns <- lahavedata$allreturns %>%
  mutate(totalreturns = totalsmall + totallarge) %>%
  select(year, smolts = smolts_med, M1 = est1SW, M2 = est2SW, RS = estRS, river_name, totalreturns)


allreturns <- bind_rows(qcreturns, nlreturns, nbreturns, nsreturns)

saveRDS(allreturns, file = "data/allreturns.rds")

