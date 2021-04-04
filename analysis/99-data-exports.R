# Data export to share with Maria Tirronen for change point analysis
library(tidyverse)
library(R.matlab)
alldata <- readRDS("data/alldata.rds")
alldatastan <- readRDS("data/alldatastan.rds")


str(alldata)

alldata$logRS <- NULL
alldata$logis_mu <- matrix(alldatastan$logis_mu, ncol = 1)
alldata$logis_sigma <- matrix(alldatastan$logis_sigma, ncol = 1)

names(alldata)[names(alldata) == "loggrilse"] <- "logSW1"

saveRDS(alldata, "data/alldata_export.rds")

writeMat("data/alldata.mat", alldata = alldata)
writeMat("data/alldatasep.mat", 
         years = alldata$years,
         logsmolts = alldata$logsmolts,
         logsmolts_cv = alldata$logsmolts_cv,
         logSW1 = alldata$logSW1,
         logSW2 = alldata$logSW2,
         returns_cv = alldata$returns_cv,
         nyears = alldata$nyears,
         river_name = alldata$river_name,
         logis_mu = alldata$logis_mu,
         logis_sigma = alldata$logis_sigma
         )

# In the end we didn't use MATLAB export and went with a simple csv table instead:

returnsdata <- alldatastan[names(alldatastan) %in% 
                             c("years", "logsmolts", "logsmolts_cv", "loggrilse", "logSW2",
                               "logSW1_cv", "logSW2_cv")] %>% do.call(cbind, .) %>% as_tibble

otherdata <- tibble(nyears = alldatastan$nyears,
       river_name = alldatastan$river_name,
       logis_mu = alldatastan$logis_mu,
       logis_sigma = alldatastan$logis_sigma,
       ny = nyears) %>% uncount(ny)

alldata_table <- bind_cols(returnsdata, otherdata) %>%
  select(river_name, year = years, logsmolts, logsmolts_cv, 
         logSW1 = loggrilse, logSW2, logSW1_cv, logSW2_cv, everything())

write.csv(alldata_table, file = "data/alldata_table.csv")


# Save csvs for sharing data with colleagues

# S1 median estimates for regime shift analysis (Anna)
s1quant <- readRDS("data/s1quant.rds")

s1medians <- s1quant %>%
  select(river_name, year, s1_median = median)

write_csv(s1medians, "data/s1medians.csv")

s1quant %>%
  select(river_name, year, s1_median = median) %>%
  mutate(z1_median = -log(s1_median))

