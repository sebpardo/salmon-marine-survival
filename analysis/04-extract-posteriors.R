library(tidyverse)
library(rstan)
library(bayesplot)

# loading stanfit objects
fithierncp <- readRDS("data/stanfit_hier_ncp.rds") # hierarchical
#fitsimplencp <- readRDS("data/stanfit_simple_ncp.rds")   # non-hierarchical

# loading data
alldatastan <- readRDS("data/alldatastan.rds")

# Indexing positions of years for each river
indexdf <- tibble(pos = as.character(1:length(alldatastan$years)),
                  river_name = rep(alldatastan$river_name, 
                                   times = alldatastan$nyears), 
                  year = alldatastan$years)

# For indexing population level posteriors
riverindexdf <- tibble(pos = as.character(1:length(alldatastan$river_name)),
                       river_name = alldatastan$river_name)

# Extracting posterior distributions of S1
s1post <- as.matrix(fithierncp, pars = "S1") %>%
  as_tibble %>%
  gather(parameter, value) %>%
  mutate(pos = str_extract(parameter, "(?<=\\[).+?(?=\\])"),
         param = str_extract(parameter, ".+?(?=\\[)")) %>%
  left_join(indexdf, by = "pos")

# These two should fail now as these years have been removed upstream
s1post %>%
  filter((river_name == "Saint-Jean River" & year == 1998))# %>%
#summarise(parameter = unique(parameter))# S1[9]
s1post %>%  
  filter((river_name == "Trinité River" & year == 2007)) # %>%
# summarise(parameter = unique(parameter)) # S1[48]

saveRDS(s1post, file = "data/s1-posteriors.rds")

s1quant <- s1post %>%
  group_by(pos, river_name, year) %>%
  summarise(median = median(value),
            q05 = quantile(value, 0.05),
            q25 = quantile(value, 0.25),
            q75 = quantile(value, 0.75),
            q95 = quantile(value, 0.95))

# Removing years with missing smolt data (which was filled with average smolt counts) 
s1quant <- s1quant %>%
  #filter(!(river_name == "Saint-Jean River" & year == 1998))  %>%
  #filter(!(river_name == "Trinité River" & year == 2007)) %>%
  group_by(river_name) %>%
  arrange(as.numeric(pos)) %>%
  mutate(zscore = scale(median))

saveRDS(s1quant, file = "data/s1quant.rds")

# Extracting posterior distributions of S2 and Pr
s2post <- as.matrix(fithierncp, pars = "S2") %>%
  as_tibble %>%
  gather(parameter, value) %>%
  mutate(pos = str_extract(parameter, "(?<=\\[).+?(?=\\])"),
         param = str_extract(parameter, ".+?(?=\\[)")) %>%
  left_join(indexdf, by = "pos")

saveRDS(s2post, file = "data/s2-posteriors.rds")

s2quant <- s2post %>%
  group_by(pos, river_name, year) %>%
  summarise(median = median(value),
            q05 = quantile(value, 0.05),
            q25 = quantile(value, 0.25),
            q75 = quantile(value, 0.75),
            q95 = quantile(value, 0.95)) %>%
  group_by(river_name) %>%
  arrange(as.numeric(pos)) %>%
  mutate(zscore = scale(median))

saveRDS(s2quant, file = "data/s2quant.rds")

prpost <- as.matrix(fithierncp, pars = "Pr") %>%
  as_tibble %>%
  gather(parameter, value) %>%
  mutate(pos = str_extract(parameter, "(?<=\\[).+?(?=\\])"),
         param = str_extract(parameter, ".+?(?=\\[)")) %>%
  left_join(indexdf, by = "pos")

saveRDS(prpost, file = "data/pr-posteriors.rds")

prquant <- prpost %>%
  group_by(pos, river_name, year) %>%
  summarise(median = median(value),
            q05 = quantile(value, 0.05),
            q25 = quantile(value, 0.25),
            q75 = quantile(value, 0.75),
            q95 = quantile(value, 0.95)) %>%
  group_by(river_name) %>%
  arrange(as.numeric(pos)) %>%
  mutate(zscore = scale(median))

saveRDS(prquant, file = "data/prquant.rds")

# Extracting logitPr estimates

logitprmupost <- as.matrix(fithierncp, pars = "logitPr_mu") %>%
  as_tibble %>%
  gather(parameter, value) %>%
  mutate(pos = str_extract(parameter, "(?<=\\[).+?(?=\\])"),
         param = str_extract(parameter, ".+?(?=\\[)")) %>%
  left_join(riverindexdf, by = "pos")

logitprsigmapost <- as.matrix(fithierncp, pars = "logitPr_sigma") %>%
  as_tibble %>%
  gather(parameter, value) %>%
  mutate(pos = str_extract(parameter, "(?<=\\[).+?(?=\\])"),
         param = str_extract(parameter, ".+?(?=\\[)")) %>%
  left_join(riverindexdf, by = "pos")

prmupost <- as.matrix(fithierncp, pars = "Pr_mu") %>%
  as_tibble %>%
  gather(parameter, value) %>%
  mutate(pos = str_extract(parameter, "(?<=\\[).+?(?=\\])"),
         param = str_extract(parameter, ".+?(?=\\[)")) %>%
  left_join(riverindexdf, by = "pos")



# Calculating CV in logitPr
logitprcvquant <- logitprmupost %>%
  rename(logitPrmu = value) %>%
  mutate(logitPrsigma = logitprsigmapost$value,
         logitPrCV = logitPrsigma/abs(logitPrmu)) %>%
  select(pos, river_name, value = logitPrCV) %>%
  mutate(param = "logitPr_CV") %>%
  group_by(param, pos, river_name) %>%
  summarise(median = median(value),
            q05 = quantile(value, 0.05),
            q25 = quantile(value, 0.25),
            q75 = quantile(value, 0.75),
            q95 = quantile(value, 0.95)) %>%
  group_by(river_name) %>%
  arrange(param, as.numeric(pos))


prmuquant <- bind_rows(logitprmupost, logitprsigmapost, prmupost) %>%
  group_by(param, pos, river_name) %>%
  summarise(median = median(value),
            q05 = quantile(value, 0.05),
            q25 = quantile(value, 0.25),
            q75 = quantile(value, 0.75),
            q95 = quantile(value, 0.95)) %>%
  group_by(river_name) %>%
  arrange(param, as.numeric(pos)) %>%
  bind_rows(logitprcvquant)


saveRDS(prmuquant, file = "data/prmuquant.rds")

saveRDS(logitprcvquant, file = "data/logitprcvquant.rds")

############################

# extracting estimates of logsmolts_true
logsmoltspost <- as.matrix(fithierncp, pars = "logsmolts_true") %>%
  as_tibble %>%
  gather(parameter, value) %>%
  mutate(pos = str_extract(parameter, "(?<=\\[).+?(?=\\])"),
         param = str_extract(parameter, ".+?(?=\\[)")) %>%
  left_join(indexdf, by = "pos")

saveRDS(logsmoltspost, file = "data/logsmoltstrue-posteriors.rds")


z1post <- as.matrix(fithierncp, pars = "Z1") %>%
  as_tibble %>%
  gather(parameter, value) %>%
  mutate(pos = str_extract(parameter, "(?<=\\[).+?(?=\\])"),
         param = str_extract(parameter, ".+?(?=\\[)")) %>%
  left_join(indexdf, by = "pos")

z2post <- as.matrix(fithierncp, pars = "Z2") %>%
  as_tibble %>%
  gather(parameter, value) %>%
  mutate(pos = str_extract(parameter, "(?<=\\[).+?(?=\\])"),
         param = str_extract(parameter, ".+?(?=\\[)")) %>%
  left_join(indexdf, by = "pos")

saveRDS(z1post, file = "data/z1-posteriors.rds")
saveRDS(z2post, file = "data/z2-posteriors.rds")

