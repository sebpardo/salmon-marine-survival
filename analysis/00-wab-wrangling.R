### WAB data
source("analysis/99-hypergeometric-bootstrap-function.R")
library(tidyverse)

wabsmolts <- read_csv("raw-data/NL-smolt-estimates-Kelly-et-al-2017-ResDoc.csv",
         na = c("-", "NA")) %>%
  rename(year = Year) %>%
  gather(river_name, smolt, -year) %>%
  filter(river_name == "Western Arm Brook")

# Loading estimated returns of 1SW, 2Sw, and RS
wab <- readRDS(file = "~/Dropbox/salmon-lh-variation/output/NL/wab-estimated-returns.rds")

# Loading raw scale age data
wabraw <- readRDS("~/Dropbox/salmon-lh-variation/output/NL/wabraw.rds") %>%
  select(year, sl, M1, M2p, RStotal)


wabwide <- pivot_wider(wabraw, names_from = sl, values_from = c(M1, M2p, RStotal)) 
wabwide[is.na(wabwide)] <- 0 # replacing NAs with zeros

# Editing this year as there were no large fish reported in the returns, yet there was
# a single fish who was scale sampled and measured 64.8 cm and was still labelled as a grilse
wabwide[wabwide$year == 1990, "M1_small"] <- 47 # Originally 46
wabwide[wabwide$year == 1990, "M1_large"] <- 0 # Originally 1
wabwide[wabwide$year == 1995, "RStotal_large"] <- 30 # originally 31, one removed to match total large return count
wabjoin <- left_join(wab, wabwide, by = "year")

wabscales <- wabjoin %>%
  mutate_all(replace_na, 0) %>%
  rename(n_small = small, 
         n_large = large, 
         scales_small_1SW = M1_small, 
         scales_small_2SW = M2p_small,
         scales_small_other = RStotal_small,
         scales_large_1SW = M1_large, 
         scales_large_2SW = M2p_large, 
         scales_large_other = RStotal_large)


wabsds <-  pmap(wabscales, hyper_boot, debug = FALSE) %>% bind_rows %>%
  select(year, starts_with("sd"))

# replacing NaNs (from years without samples) with mean value
wabsds[is.nan(wabsds$sdlog2SW),"sdlog2SW"] <- 0.99 #  mean(wabsds$sdlog2SW, na.rm = TRUE)


wabepsilon <- select(wabscales, year, n_small, n_large, totalreturns, scales_small_1SW:scales_large_other) %>%
left_join(wabsds, by = "year") %>%
  mutate(river_name = "Western Arm Brook") %>%
  select(river_name, everything())

saveRDS(wabepsilon, file = "data/wabepsilon.rds")

waball <- right_join(wab, wabsmolts, by = "year") %>%
  left_join(wabsds, by = "year") %>%
  select(year, smolt, estM1, estM2p, estRS, totalreturns, starts_with("sd"))


wabreturns <- waball %>%
  mutate(smoltlag = lag(smolt),
         #smolt025lag = lag(smolt_025),
         #smolt975lag = lag(smolt_975),
         est_1SW = estM1,
         est_2SWlead = lead(estM2p),
         est_RSlead = lead(estRS),
         sdlog2SWlead = lead(sdlog2SW)
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
                logSW1_cv = wabreturns$sdlog1SW,
                logSW2_cv = wabreturns$sdlog2SWlead,
                  nyears = length(wabreturns$year),
                  river_name = "Western Arm Brook",
                  allreturns = waball %>% mutate(river_name = "Western Arm Brook")
) 

saveRDS(wabdata, file = "data/wabdata.rds")

