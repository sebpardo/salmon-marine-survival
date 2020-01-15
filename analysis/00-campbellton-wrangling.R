### WAB data

library(tidyverse)

campsmolts <- read_csv("raw-data/NL-smolt-estimates-Kelly-et-al-2017-ResDoc.csv",
         na = c("-", "NA")) %>%
  rename(year = Year) %>%
  gather(river_name, smolt, -year) %>%
  filter(river_name == "Campbellton River") %>%
  na.omit()

# Loading estimated returns of 1SW, 2Sw, and RS
camp <- readRDS(file = "~/Dropbox/salmon-lh-variation/output/NL/campbellton-estimated-returns.rds")

campall <- right_join(camp, campsmolts, by = "year") %>%
  select(year, smolt, estM1, estM2p, estRS, totalreturns) %>%
  mutate(river_name = "Campbellton River")
  


campreturns <- campall %>%
  mutate(smoltlag = lag(smolt),
         #smolt025lag = lag(smolt_025),
         #smolt975lag = lag(smolt_975),
         est_1SW = estM1,
         est_2SWlead = lead(estM2p),
         est_RSlead = lead(estRS)
         #logsmoltsdlag = sqrt(300) * (log(smolt975lag) - log(smolt025lag))/3.92
  ) %>%
  filter(year %in% 1994:2015)

campreturns %>%
  filter(est_2SWlead != 0) %>%
  summarise(min2SW = min(est_2SWlead)/2)

ggplot(campreturns, aes(year, smolt)) + geom_point() + geom_line() # +

campreturns[campreturns$est_2SWlead == 0, "est_2SWlead"] <- 1.351 # roughly min(est_2SW)/2

campreturn_rate <- campreturns %>%
  mutate(return_rate = est_1SW/smoltlag) %>%
  select(year, return_rate) %>%
  mutate(river_name = "Campbellton River")

ggplot(campreturn_rate, aes(year, return_rate)) +
  geom_point() + geom_line() + ggtitle("Campbellton return rates")

saveRDS(campreturn_rate, file = "data/campbelltonreturnrate.rds")

# Creating list for use in TMB/Stan
campdata <- list(years = campreturns$year,
                  logsmolts = log(campreturns$smoltlag),
                  #logsmolts_SE = nashwaakreturns$logsmoltsdlag,
                  logsmolts_cv =  rep(0.05, length(campreturns$smoltlag)), # should be replaced with actual uncertainty, if available
                  loggrilse = log(campreturns$est_1SW),
                  logSW2 = log(campreturns$est_2SWlead),
                  logRS = log(campreturns$est_RSlead),
                  #N =  length(wabreturns$small[-1]),
                  #Pr = rep(0.8, nrow(wabreturns)))
                  returns_cv = 0.01, # ASSUMPTION, CV of RETURN abundance counts
                  nyears = length(campreturns$year),
                  river_name = "Campbellton River",
                  allreturns = campall
) 

saveRDS(campdata, file = "data/campbelltondata.rds")
