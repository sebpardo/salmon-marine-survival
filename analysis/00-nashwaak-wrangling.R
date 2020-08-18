# Nashwaak data

library(tidyverse)

# Load Naskwaak data from salmon-inverse-matrix repo
load("~/Dropbox/salmon-inverse-matrix/olddata/nashwaak.rda")

# We use the return data frame that includes wild salmon only, as it
nashwaak_wild_returns
nashwaak_returns_wild_hatch
nashwaak_smolts
# Smolt data from 1998 to 2016

# Saving filtered data frame in new object
nashwaakdat <- right_join(
  select(nashwaak_smolts, year, smolt_mode, smolt_025, smolt_975),
  select(nashwaak_wild_returns, year, est_1SW, est_2SW, est_RS), 
  by = "year")


# Missing data:
#
# Smolt data missing for 2017 and 2018 
# so after lag in smolts just need to exclude year 2018

nashwaakreturns <- nashwaakdat %>%
  mutate(smoltdiff =  ((smolt_mode-smolt_025)+(smolt_975-smolt_mode))/2,
         smoltdiff_upper =  smolt_975-smolt_mode,
         smoltdiff_lower =  smolt_mode-smolt_025,
         smoltsd = smoltdiff/1.96,
         smoltscv = smoltsd/smolt_mode,
         smoltlag = lag(smolt_mode),
         smolt025lag = lag(smolt_025),
         smolt975lag = lag(smolt_975),
         est_2SWlead = lead(est_2SW),
         est_RSlead = lead(est_RS),
         smoltscvlag = lag(smoltscv)) %>% 
  filter(year %in% 1999:2017)


nashwaakreturns

ggplot(nashwaakreturns, aes(year, smolt_mode)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin = smolt_025, ymax = smolt_975))

# Creating list for use in TMB
nashwaakdata <- list(years = nashwaakreturns$year,
  logsmolts = log(nashwaakreturns$smoltlag),
  logsmolts_cv = nashwaakreturns$smoltscvlag,
  loggrilse = log(nashwaakreturns$est_1SW),
  logSW2 = log(nashwaakreturns$est_2SWlead),
  logRS = log(nashwaakreturns$est_RSlead),
  #N =  length(connereturns$small[-1]),
  #Pr = rep(0.8, nrow(connereturns)))
  returns_cv = 0.01, # ASSUMPTION, CV of return abundance counts
  logSW1_cv = rep(0.01, length(nashwaakreturns$year)), # NO DATA SO USING VALUES BASED ON MEAN TRINITE ESTIMATES
  logSW2_cv =  rep(0.06, length(nashwaakreturns$year)),
  nyears = length(nashwaakreturns$year),
  river_name = "Nashwaak River",
  allreturns = nashwaakdat  %>% filter(year >= 1980) %>% mutate(river_name = "Nashwaak River")
) 

saveRDS(nashwaakdata, file = "data/nashwaakdata.rds")
