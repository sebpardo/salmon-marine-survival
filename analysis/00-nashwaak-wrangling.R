# Load Naskwaak data from salmon-inverse-matrix repo
load("~/Dropbox/salmon-inverse-matrix/data/nashwaak.rda")

# We use the return data frame that includes wild salmon only, as it
nashwaak_returns_wild_hatch

# Saving filtered data frame in new object
nashwaakdat <- right_join(
  select(nashwaak_smolts, year, smolt_mode, smolt_025, smolt_975),
  select(nashwaak_wild_returns, year, est_1SW, est_2SW, est_RS), 
  by = "year")



nashwaakreturns <- nashwaakdat %>%
  mutate(smoltlag = lag(smolt_mode),
         smolt025lag = lag(smolt_025),
         smolt975lag = lag(smolt_975),
         est_2SWlead = lead(est_2SW),
         est_RSlead = lead(est_RS),
         logsmoltsdlag = sqrt(300) * (log(smolt975lag) - log(smolt025lag))/3.92) %>% 
  filter(year %in% 1999:2017)

ggplot(nashwaakreturns, aes(year, smolt_mode)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin = smolt_025, ymax = smolt_975))

# Creating list for use in TMB
nashwaakdata <- list(years = nashwaakreturns$year,
  logsmolts = log(nashwaakreturns$smoltlag),
  logsmolts_SE = nashwaakreturns$logsmoltsdlag,
  #logsmolts_SE =  rep(0.05, length(nashwaakreturns$smoltlag)), # should be replaced with actual uncertainty
  loggrilse = log(nashwaakreturns$est_1SW),
  logSW2 = log(nashwaakreturns$est_2SWlead),
  logRS = log(nashwaakreturns$est_RSlead),
  #N =  length(connereturns$small[-1]),
  #Pr = rep(0.8, nrow(connereturns)))
  cv = 0.01 # ASSUMPTION, CV of return abundance counts
) 
