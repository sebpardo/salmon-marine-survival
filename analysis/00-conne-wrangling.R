# Conne data

library(tidyverse)

# Loading estimated returns of 1SW, 2Sw, and RS
conne <- readRDS(file = "~/Dropbox/salmon-lh-variation/output/NL/conne-estimated-returns.rds")

# Loading small, large, and smolt data
connewithsmolts <- readRDS(file = "~/Dropbox/salmon-inverse-matrix/data/connereturns.rds")

# joining by year
conneall <- left_join(conne, connewithsmolts, by = "year") %>% print(n = 32)
  select(year, smolt, estM1, estM2p, estRS, totalreturns)

connereturns <- conneall %>%
  mutate(smoltlag = lag(smolt),
         #smolt025lag = lag(smolt_025),
         #smolt975lag = lag(smolt_975),
         est_1SW = estM1,
         est_2SWlead = lead(estM2p),
         est_RSlead = lead(estRS),
         #logsmoltsdlag = sqrt(300) * (log(smolt975lag) - log(smolt025lag))/3.92
         ) %>% 
 filter(year %in% 1999:2015) #  2017) # NEED 2017 SCALE DATA FOR 1SW/2SW/RS ESTIMATION IN BOTH 2016 AND 2017, ASK BRIAN

ggplot(connereturns, aes(year, smolt)) + geom_point() + geom_line() # +
 # geom_errorbar(aes(ymin = smolt_025, ymax = smolt_975))

# Removing zeros in 2SW-RS estimates as they cannot be log transformed
connereturns$est_2SWlead[connereturns$est_2SWlead == 0] <- 0.0001

print(connereturns, n = 32)

# Creating list for use in TMB
connedata <- list(years = connereturns$year,
                     logsmolts = log(connereturns$smoltlag),
                     #logsmolts_SE = nashwaakreturns$logsmoltsdlag,
                     logsmolts_SE =  rep(1, length(connereturns$smoltlag)), # should be replaced with actual uncertainty
                     loggrilse = log(connereturns$est_1SW),
                     logSW2 = log(connereturns$est_2SWlead),
                     logRS = log(connereturns$est_RSlead),
                     #N =  length(connereturns$small[-1]),
                     #Pr = rep(0.8, nrow(connereturns)))
                     cv = 0.01, # ASSUMPTION, CV of return abundance counts
                  nyears = length(connereturns$year),
                  river_name = "Conne River"
                  ) 

saveRDS(connedata, file = "data/connedata.rds")
