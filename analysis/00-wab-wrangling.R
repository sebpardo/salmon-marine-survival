### WAB data

library(tidyverse)

wabsmolts <- read_csv("raw-data/NL-smolt-estimates-Kelly-et-al-2017-ResDoc.csv",
         na = c("-", "NA")) %>%
  rename(year = Year) %>%
  gather(river_name, smolt, -year) %>%
  filter(river_name == "Western Arm Brook")

# Loading estimated returns of 1SW, 2Sw, and RS
wab <- readRDS(file = "~/Dropbox/salmon-lh-variation/output/NL/wab-estimated-returns.rds")

waball <- right_join(wab, wabsmolts, by = "year") %>%
  select(year, smolt, estM1, estM2p, estRS, totalreturns)


wabreturns <- waball %>%
  mutate(smoltlag = lag(smolt),
         #smolt025lag = lag(smolt_025),
         #smolt975lag = lag(smolt_975),
         est_1SW = estM1,
         est_2SWlead = lead(estM2p),
         est_RSlead = lead(estRS)
         #logsmoltsdlag = sqrt(300) * (log(smolt975lag) - log(smolt025lag))/3.92
  ) %>%
  filter(year %in% 1974:2015)

wabreturns %>%
  filter(est_2SWlead != 0) %>%
  summarise(min2SW = min(est_2SWlead))

ggplot(wabreturns, aes(year, smolt)) + geom_point() + geom_line() # +

wabreturns[wabreturns$est_2SWlead == 0, "est_2SWlead"] <- 0.1 # roughly min(est_2SW)/2

wabreturn_rate <- wabreturns %>%
  mutate(return_rate = est_1SW/smoltlag) %>%
  select(year, return_rate) %>%
  mutate(river_name = "Western Arm Brook")

ggplot(wabreturn_rate, aes(year, return_rate)) + 
  geom_point() + geom_line() + ggtitle("WAB return rates")


saveRDS(wabreturn_rate, file = "data/wabreturnrate.rds")

# Creating list for use in TMB
wabdata <- list(years = wabreturns$year,
                  logsmolts = log(wabreturns$smoltlag),
                  #logsmolts_SE = nashwaakreturns$logsmoltsdlag,
                  logsmolts_cv =  rep(0.05, length(wabreturns$smoltlag)), # should be replaced with actual uncertainty, if available
                  loggrilse = log(wabreturns$est_1SW),
                  logSW2 = log(wabreturns$est_2SWlead),
                  logRS = log(wabreturns$est_RSlead),
                  #N =  length(wabreturns$small[-1]),
                  #Pr = rep(0.8, nrow(wabreturns)))
                  returns_cv = 0.01, # ASSUMPTION, CV of RETURN abundance counts
                  nyears = length(wabreturns$year),
                  river_name = "Western Arm Brook",
                  allreturns = waball %>% mutate(river_name = "Western Arm Brook")
) 

saveRDS(wabdata, file = "data/wabdata.rds")

