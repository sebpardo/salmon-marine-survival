### Saint-Jean Data

library(tidyverse)

saintjean <- readRDS(file = "~/Dropbox/salmon-bayesian-survival/qc-stjean-trinite-all.rds") %>%
  filter(river == "Saint-Jean", year >= 1981) %>%
  spread(stage, numbers) %>%
  select(year, smolts, returns_grilse, returns_2sw_est, returns_3sw_est, returns_rs_est)# %>%
  #filter(year %in% 1998:2014)


saintjeanreturns <- saintjean %>%
  mutate(smoltlag = lag(smolts),
         #smolt025lag = lag(smolt_025),
         #smolt975lag = lag(smolt_975),
         est_1SW = returns_grilse,
         est_2SWlead = lead(returns_2sw_est),
         est_RSlead = lead(returns_rs_est),
         #logsmoltsdlag = sqrt(300) * (log(smolt975lag) - log(smolt025lag))/3.92
  ) %>% 
  filter(year %in% 1990:2014)


# Removing zeros in 2SW-RS estimates as they cannot be log transformed
saintjeanreturns$est_RSlead[saintjeanreturns$est_RSlead == 0] <- 0.0001


# using average of contiguous 4 years (+-2) instead
saintjeanreturns[is.na(saintjeanreturns$smoltlag), "smoltlag"] <-  mean(saintjeanreturns$smoltlag[saintjeanreturns$year %in% 1996:2000], na.rm = TRUE)


ggplot(saintjeanreturns, aes(year, smolts)) + geom_point() + geom_line() # +
# geom_errorbar(aes(ymin = smolt_025, ymax = smolt_975))

# Creating list for use in TMB
saintjeandata <- list(years = saintjeanreturns$year,
                    logsmolts = log(saintjeanreturns$smoltlag),
                    #logsmolts_SE = nashwaakreturns$logsmoltsdlag,
                    logsmolts_SE =  rep(1, length(saintjeanreturns$smoltlag)), # should be replaced with actual uncertainty
                    loggrilse = log(saintjeanreturns$est_1SW),
                    logSW2 = log(saintjeanreturns$est_2SWlead),
                    logRS = log(saintjeanreturns$est_RSlead),
                    #N =  length(saintjeanreturns$small[-1]),
                    #Pr = rep(0.8, nrow(connereturns)))
                    cv = 0.01, # ASSUMPTION, CV of return abundance counts
                    nyears = length(saintjeanreturns$year),
                    river_name = "Saint-Jean River"
) 


saveRDS(saintjeandata, file = "data/saintjeandata.rds")
