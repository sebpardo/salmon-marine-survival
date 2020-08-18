# Conne data

source("analysis/99-hypergeometric-bootstrap-function.R")
library(tidyverse)


# Loading estimated returns of 1SW, 2Sw, and RS
conne <- readRDS(file = "~/Dropbox/salmon-lh-variation/output/NL/conne-estimated-returns.rds")

# Loading raw scale age data
conneraw <- readRDS("~/Dropbox/salmon-lh-variation/output/NL/conneraw.rds") %>%
  select(year, sl, M1, M2p, RStotal)

# Loading small, large, and smolt data
# Smolt years missing: up to 1988, 2016 and 2017, which
# is the duration of the dataset, so no need to remove years in the middle
connewithsmolts <- readRDS(file = "~/Dropbox/salmon-inverse-matrix/olddata/connereturns.rds")

# Adding missing years based on Dempson et al 2003 (Ch 1)
connewithsmolts[connewithsmolts$year == 1987, "smolt"] <- 74585
connewithsmolts[connewithsmolts$year == 1988, "smolt"] <- 65692


connewide <- pivot_wider(conneraw, names_from = sl, values_from = c(M1, M2p, RStotal)) 
connewide[is.na(connewide)] <- 0 # replacing NAs with zeros

connejoin <- left_join(connewithsmolts, connewide, by = "year") %>%
   filter(year >= 1986 & year < 2017) %>% 
  select(-smolt) 

connescales <- connejoin %>%
  rename(n_small = small, 
         n_large = large, 
         scales_small_1SW = M1_small, 
         scales_small_2SW = M2p_small,
         scales_small_other = RStotal_small,
         scales_large_1SW = M1_large, 
         scales_large_2SW = M2p_large, 
         scales_large_other = RStotal_large)


connesds <-  pmap(connescales, hyper_boot) %>% bind_rows %>%
  select(year, starts_with("sd"))

# replacing NaNs (from years without samples) with mean value
connesds[is.nan(connesds$sdlog2SW),"sdlog2SW"] <- 0.99 #  mean(wabsds$sdlog2SW, na.rm = TRUE)

conneepsilon <- connescales %>%
  mutate(totalreturns = n_small + n_large) %>%
  select(year, n_small, n_large, totalreturns, scales_small_1SW:scales_large_other) %>%
  left_join(connesds, by = "year") %>%
  mutate(river_name = "Conne River") %>%
  select(river_name, everything())

saveRDS(conneepsilon, file = "data/conneepsilon.rds")



# joining by year
conneall <- left_join(conne, connewithsmolts, by = "year") %>%
  left_join(connesds, by = "year") %>%
  select(year, smolt, estM1, estM2p, estRS, totalreturns, starts_with("sd")) %>%
  mutate(river_name = "Conne River")
  

connereturns <- conneall %>%
  mutate(smoltlag = lag(smolt),
         #smolt025lag = lag(smolt_025),
         #smolt975lag = lag(smolt_975),
         est_1SW = estM1,
         est_2SWlead = lead(estM2p),
         est_RSlead = lead(estRS),
         sdlog2SWlead = lead(sdlog2SW),
         #logsmoltsdlag = sqrt(300) * (log(smolt975lag) - log(smolt025lag))/3.92
         ) %>% 
 filter(year %in% 1988:2015) #  2017) # NEED 2017 SCALE DATA FOR 1SW/2SW/RS ESTIMATION IN BOTH 2016 AND 2017, ASK BRIAN

ggplot(connereturns, aes(year, smolt)) + geom_point() + geom_line() # +
 # geom_errorbar(aes(ymin = smolt_025, ymax = smolt_975))

# Removing zeros in 2SW-RS estimates as they cannot be log transformed
# not needed after fixing years with very few large scales (and 0s in estimated 2SW returns)
# connereturns$est_2SWlead[connereturns$est_2SWlead == 0] <- 0.0001

filter(connereturns, est_2SWlead != 0) %>% pull(est_2SWlead) %>% min(.)/2 # picking min(RS)/2

connereturns[connereturns$est_2SWlead == 0, "est_2SWlead"] <- 3.7 # roughly min(est_2SW)/2



print(connereturns, n = 32)

connereturn_rate <- connereturns %>%
  mutate(return_rate = est_1SW/smoltlag) %>%
  select(year, return_rate) %>%
  mutate(river_name = "Conne River")

ggplot(connereturn_rate, aes(year, return_rate)) +
  geom_point() + geom_line() + ggtitle("Conne return rates")

saveRDS(connereturn_rate, file = "data/connereturnrate.rds")

# Creating list for use in TMB
connedata <- list(years = connereturns$year,
                     logsmolts = log(connereturns$smoltlag),
                     #logsmolts_SE = nashwaakreturns$logsmoltsdlag,
                     logsmolts_cv =  rep(0.1, length(connereturns$smoltlag)), # should be replaced with actual uncertainty but CV of log data ~= SD
                     loggrilse = log(connereturns$est_1SW),
                     logSW2 = log(connereturns$est_2SWlead),
                     logRS = log(connereturns$est_RSlead),
                     #N =  length(connereturns$small[-1]),
                     #Pr = rep(0.8, nrow(connereturns)))
                     returns_cv = 0.01, # ASSUMPTION, CV of return abundance counts
                     logSW1_cv = connereturns$sdlog1SW,
                  logSW2_cv = connereturns$sdlog2SWlead,
                  nyears = length(connereturns$year),
                  river_name = "Conne River",
                  allreturns = conneall
                  ) 

saveRDS(connedata, file = "data/connedata.rds")
