### Trinite Data

library(tidyverse)

trinite <- readRDS(file = "~/Dropbox/salmon-bayesian-survival/qc-stjean-trinite-all.rds") %>%
  filter(river == "Trinité") %>%
  spread(stage, numbers) %>%
  select(year, smolts, returns_grilse, returns_2sw_est, returns_3sw_est, returns_rs_est) # %>%
 # filter(year %in% 1984:2017)

# NEED TO ADD YEARS AFTER 2006

trinite %>%
  mutate(smolts = smolts/50) %>%
  gather(stage, number, - year) %>%
  ggplot(aes(year, number, color = stage)) + geom_point() +
  geom_line() +
  scale_y_continuous("return number", sec.axis = sec_axis(~ . * 50, name = "smolt number"))


trinitereturns <- trinite %>%
  mutate(smoltlag = lag(smolts),
         #smolt025lag = lag(smolt_025),
         #smolt975lag = lag(smolt_975),
         est_1SW = returns_grilse,
         est_2SWlead = lead(returns_2sw_est),
         est_RSlead = lead(returns_rs_est),
         #logsmoltsdlag = sqrt(300) * (log(smolt975lag) - log(smolt025lag))/3.92
        ) %>% 
  filter(year %in% 1985:2016) %>% # lead 2017 is missing 2SW and RS
  filter(year != 2007) # Removing smolt year 2006 as its missing smolt estimates 

# YEAR WITH NO SMOLT DATA (NOT NEEDED ANYMORE AS AM EXCLUDING THOSE YEARS)
# trinitedata$smolts[trinitedata$year == 2007] <- 1
# Can't use a dummy value since it will screw up with resulting 
# distributions for the population
# using average of contiguous 4 years (+-2) instead
# trinitereturns[is.na(trinitereturns$smoltlag), "smoltlag"] <-  mean(trinitereturns$smoltlag[trinitereturns$year %in% 2005:2009], na.rm = TRUE)




ggplot(trinitereturns, aes(year, smolts)) + geom_point() + geom_line() +
    geom_point(aes(year, smoltlag), color = "red") # +
# geom_errorbar(aes(ymin = smolt_025, ymax = smolt_975))

# Creating list for use in TMB
trinitedata <- list(years = trinitereturns$year,
                  logsmolts = log(trinitereturns$smoltlag),
                  #logsmolts_SE = nashwaakreturns$logsmoltsdlag,
                  logsmolts_cv =  rep(0.1, length(trinitereturns$smoltlag)), # should be replaced with actual uncertainty
                  loggrilse = log(trinitereturns$est_1SW),
                  logSW2 = log(trinitereturns$est_2SWlead),
                  logRS = log(trinitereturns$est_RSlead),
                  #N =  length(connereturns$small[-1]),
                  #Pr = rep(0.8, nrow(connereturns)))
                  returns_cv = 0.01, # ASSUMPTION, CV of return abundance counts
                  nyears = length(trinitereturns$year),
                  river_name = "Trinité River",
                  allreturns = trinite %>% mutate(river_name = "Trinité River")
                  ) 

saveRDS(trinitedata, file = "data/trinitedata.rds")
