### WAB data

source("analysis/99-hypergeometric-bootstrap-function.R")
library(tidyverse)

campsmolts <- read_csv("raw-data/NL-smolt-estimates-Kelly-et-al-2017-ResDoc.csv",
         na = c("-", "NA")) %>%
  rename(year = Year) %>%
  gather(river_name, smolt, -year) %>%
  filter(river_name == "Campbellton River") %>%
  na.omit()

# Loading estimated returns of 1SW, 2Sw, and RS
camp <- readRDS(file = "~/Dropbox/salmon-lh-variation/output/NL/campbellton-estimated-returns.rds")

# Loading raw scale age data
campraw <- readRDS("~/Dropbox/salmon-lh-variation/output/NL/campraw.rds") %>%
  select(year, sl, M1, M2p, RStotal)

campwide <- pivot_wider(campraw, names_from = sl, values_from = c(M1, M2p, RStotal)) 
campwide[is.na(campwide)] <- 0 # replacing NAs with zeros

campjoin <- left_join(camp, campwide, by = "year")

campscales <- campjoin %>%
  mutate_all(replace_na, 0) %>%
  rename(n_small = small, 
         n_large = large, 
         scales_small_1SW = M1_small, 
         scales_small_2SW = M2p_small,
         scales_small_other = RStotal_small,
         scales_large_1SW = M1_large, 
         scales_large_2SW = M2p_large, 
         scales_large_other = RStotal_large)


campsds <-  pmap(campscales, hyper_boot) %>% bind_rows %>%
  select(year, starts_with("sd"))

# replacing NaNs (from years without samples) with mean value
campsds[is.nan(campsds$sdlog2SW),"sdlog2SW"] <- 0.99 #  mean(wabsds$sdlog2SW, na.rm = TRUE)


campepsilon <- select(campscales, year, n_small, n_large, totalreturns, scales_small_1SW:scales_large_other) %>%
  left_join(campsds, by = "year") %>%
  mutate(river_name = "Campbellton River") %>%
  select(river_name, everything())

saveRDS(campepsilon, file = "data/campepsilon.rds")

campscales2011 <- filter(campscales, year == 2011)
campscales2014 <- filter(campscales, year == 2014)


debug(hyper_boot)
campsds2011 <- hyper_boot(year = campscales2011$year,
         campscales2011$n_small,
         campscales2011$n_large,
         campscales2011$scales_small_1SW,
         campscales2011$scales_small_2SW,
         campscales2011$scales_small_other,
         campscales2011$scales_large_1SW,
         campscales2011$scales_large_2SW,
         campscales2011$scales_large_other,
         iter = 10000, list = TRUE)

campsds2014 <- hyper_boot(year = campscales2014$year,
                          campscales2014$n_small,
                          campscales2014$n_large,
                          campscales2014$scales_small_1SW,
                          campscales2014$scales_small_2SW,
                          campscales2014$scales_small_other,
                          campscales2014$scales_large_1SW,
                          campscales2014$scales_large_2SW,
                          campscales2014$scales_large_other,
                          iter = 10000, list = TRUE)

oldpar <- par()

png("figures/camp-scale-epsilon-example.png", width = 8, height = 6, units = "in", res = 300)
par(mfrow = c(2,2), mar = c(5,4,2,2), oma = c(1,1,3,0))
hist(campsds2011$sm1 + campsds2011$la1, breaks = seq(0,5500, by = 100), col = "lightgray", 
     probability = TRUE, main = "", xlab = "Number")
mtext("a)", side = 3, adj =0.05, cex = 1.3, line = -1)
hist(log(campsds2011$sm1 + campsds2011$la1), breaks = seq(0,10, by = 0.05), col = "lightgray", 
     probability = TRUE, main = "", #Draws of small and large 1SW return abundances in log-space",
     xlab = "log Number")
mtext("b)", side = 3, adj =0.05, cex = 1.3, line = -1)
mtext(expression(paste(italic(epsilon[1]),"= 0.711")), side = 1, adj =0.1, cex = 1.2, line = -5)
#abline(v = campsds2011$meanlog1SW)
curve(dnorm(x, mean=campsds2011$meanlog1SW, sd=campsds2011$sdlog1SW), add=TRUE, n = 10001, col= "blue", lwd = 1) 

hist(campsds2014$sm1 + campsds2014$la1, breaks = seq(0,5500, by = 100), col = "lightgray", 
     probability = TRUE, main = "", xlab = "Number")
mtext("c)", side = 3, adj =0.05, cex = 1.3, line = -1)
hist(log(campsds2014$sm1 + campsds2014$la1), breaks = seq(0,10, by = 0.05), col = "lightgray", 
     probability = TRUE, main = "", #Draws of small and large 1SW return abundances in log-space",
     xlab = "log Number")
mtext("d)", side = 3, adj =0.05, cex = 1.3, line = -1)
mtext(expression(paste(italic(epsilon[1]),"= 0.033")), side = 1, adj =0.1, cex = 1.2, line = -5)
curve(dnorm(x, mean=campsds2014$meanlog1SW, sd=campsds2014$sdlog1SW), n = 10001, add=TRUE, col= "blue", lwd = 1) 

mtext("Draws of small and large 1SW return abundances in linear and log-space",outer = TRUE, 
      side = 3, cex = 1.3, line = 1)
mtext("Campbellton, year 2011, 0 1SW scale samples, 5443 total returns" ,outer = TRUE, 
      side = 3, cex = 1.1, line = -1.1, adj = 0.1)
mtext("Campbellton, year 2014, 129 1SW scale samples, 4533 total returns" ,outer = TRUE, 
      side = 1, cex = 1.1, line = -16.4, adj = 0.1)
dev.off()

h <- hist(log(campsds2011$sm1 + campsds2011$la1), breaks = seq(0,10, by = 0.5), freq=FALSE,
          col = "lightgray", xlab = "Accuracy", main = "Overall") 
xfit <- seq(min(campsds2011$sm1 + campsds2011$la1), max(campsds2011$sm1 + campsds2011$la1), length = 40) 
yfit <- dnorm(xfit, mean = mean(campsds2011$sm1 + campsds2011$la1), sd = sd(campsds2011$sm1 + campsds2011$la1)) 
yfit <- yfit * diff(h$mids[1:2]) * length(campsds2011$sm1 + campsds2011$la1) 
lines(xfit, yfit, col = "black", lwd = 2)

df2001 <- bind_rows(tibble(est = campsds2011$sm1, agestage = "sm1"),
tibble(est = campsds2011$la1,agestage = "la1")) %>%
  mutate(sum = sm1 + sm2, logsum = log(sum))

ggplot(df2001, aes(x = est, col = agestage)) + geom_histogram()


campall <- right_join(camp, campsmolts, by = "year") %>%
  left_join(campsds, by = "year") %>%
  select(year, smolt, estM1, estM2p, estRS, totalreturns, starts_with("sd")) %>%
  mutate(river_name = "Campbellton River")
  

campreturns <- campall %>%
  mutate(smoltlag = lag(smolt),
         #smolt025lag = lag(smolt_025),
         #smolt975lag = lag(smolt_975),
         est_1SW = estM1,
         est_2SWlead = lead(estM2p),
         est_RSlead = lead(estRS),
         sdlog2SWlead = lead(sdlog2SW)
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
                 logSW1_cv = campreturns$sdlog1SW,
                 logSW2_cv = campreturns$sdlog2SWlead,
                  nyears = length(campreturns$year),
                  river_name = "Campbellton River",
                  allreturns = campall
) 

saveRDS(campdata, file = "data/campbelltondata.rds")
