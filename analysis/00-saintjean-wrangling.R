### Saint-Jean Data

library(tidyverse)
library(readxl)


saintjean <- readRDS(file = "~/Dropbox/salmon-bayesian-survival/qc-stjean-trinite-all.rds") %>%
  filter(river == "Saint-Jean", year >= 1981) %>%
  spread(stage, numbers) %>%
  select(year, smolts, returns_grilse, returns_2sw_est, returns_3sw_est, returns_rs_est)# %>%
  #filter(year %in% 1998:2014)

# smolt data from 1989 onwards, missing 1997 and 2015 


sjsmolts <- read_excel("raw-data/Estimation dévalaison St-Jean-Trinité 1984-2018.xlsx", sheet = "Saint-Jean") %>%
  select("Year", "N Min (95%)":"Smolt migrate to sea") %>%
  rename(year = Year,
         n_min = `N Min (95%)`,
         n = `N Estimation`,
         n_max = `N max (95%)`,
         smolt_prod = `Smolt produce`,
         smolt_tosea = `Smolt migrate to sea`) %>%
  mutate(n_sea_diff = n - smolt_tosea,
         smolt_est = smolt_tosea,
         smolt_min = n_min - n_sea_diff,
         smolt_max = n_max - n_sea_diff)


excel_sheets("raw-data/Estimation dévalaison St-Jean-Trinité 1984-2018.xlsx")


saintjeanreturns <- saintjean %>%
  mutate(smoltlag = lag(smolts),
         #smolt025lag = lag(smolt_025),
         #smolt975lag = lag(smolt_975),
         est_1SW = returns_grilse,
         est_2SWlead = lead(returns_2sw_est),
         est_RSlead = lead(returns_rs_est),
         #logsmoltsdlag = sqrt(300) * (log(smolt975lag) - log(smolt025lag))/3.92
  ) %>% 
  filter(year %in% 1990:2015) %>% # updated from 2014 to 2015
  filter(year != 1998) # removing smotl year 1997 as its missing smolt data

# Removing zeros in 2SW-RS estimates as they cannot be log transformed
filter(saintjeanreturns, est_RSlead != 0) %>% pull(est_RSlead) %>% min(.)/2 # picking min(RS)/2
saintjeanreturns$est_RSlead[saintjeanreturns$est_RSlead == 0] <- 4.678571

# using average of contiguous 4 years (+-2) instead to fill in missing 1997 smolt count
# saintjeanreturns[is.na(saintjeanreturns$smoltlag), "smoltlag"] <-  mean(saintjeanreturns$smoltlag[saintjeanreturns$year %in% 1996:2000], na.rm = TRUE)


ggplot(saintjean, aes(year, smolts)) + geom_point() + geom_line() # +
# geom_errorbar(aes(ymin = smolt_025, ymax = smolt_975))

# Creating list for use in TMB
saintjeandata <- list(years = saintjeanreturns$year,
                    logsmolts = log(saintjeanreturns$smoltlag),
                    #logsmolts_SE = nashwaakreturns$logsmoltsdlag,
                    logsmolts_cv =  rep(0.1, length(saintjeanreturns$smoltlag)), # should be replaced with actual uncertainty
                    loggrilse = log(saintjeanreturns$est_1SW),
                    logSW2 = log(saintjeanreturns$est_2SWlead),
                    logRS = log(saintjeanreturns$est_RSlead),
                    #N =  length(saintjeanreturns$small[-1]),
                    #Pr = rep(0.8, nrow(connereturns)))
                    returns_cv = 0.01, # ASSUMPTION, CV of return abundance counts
                    nyears = length(saintjeanreturns$year),
                    river_name = "Saint-Jean River",
                    allreturns = saintjean %>% mutate(river_name = "Saint-Jean River")
) 


saveRDS(saintjeandata, file = "data/saintjeandata.rds")
