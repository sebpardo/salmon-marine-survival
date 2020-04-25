# Exploring correlations among trends of S1, smolt, and return abundances

library(tidyverse)
library(corrplot)
library(cowplot)
source("analysis/99-helper-functions.R")

fithierncp <- readRDS("data/stanfit_hier_ncp.rds") # hierarchical

alldatastan <- readRDS("data/alldatastan.rds")
s1quant <- readRDS("data/s1quant.rds")
allreturns <- readRDS("data/allreturns.rds") %>% 
  group_by(river_name) %>%
  mutate(zscoresmolts = scale(log(smolts)), zscoretotal = scale(log(totalreturns))) %>%
  ungroup()
  
estreturns <- tibble(river_name = rep(alldatastan$river_name, times = alldatastan$nyears),
       year = alldatastan$years,
       #logsmolts = alldatastan$logsmolts, 
       log1SW = alldatastan$loggrilse,
       log2SW = alldatastan$logSW2)

river_names <- rep(alldatastan$river_name, times = alldatastan$nyears)

propgrilse <- estreturns %>%
  mutate(propgrilse = exp(log1SW)/(exp(log1SW)+exp(log2SW))) %>%
  group_by(river_name) %>%
  summarise(propgrilse = mean(propgrilse)) %>%
  arrange(desc(propgrilse))

s1post <- readRDS(file = "data/s1-posteriors.rds")
s2post <- readRDS(file = "data/s2-posteriors.rds")
z1post <- readRDS(file = "data/z1-posteriors.rds")
z2post <- readRDS(file = "data/z2-posteriors.rds")
prpost <- readRDS(file = "data/pr-posteriors.rds")
logsmoltspost <- readRDS(file = "data/logsmoltstrue-posteriors.rds")

# estimating mu1 and mu2 from posteriors
mu1 <- logsmoltspost$value + log(prpost$value) - z1post$value
mu2 <- logsmoltspost$value - z1post$value + log(1 - prpost$value)  - z2post$value

z1mat <- as.matrix(fithierncp, pars = "Z1") 
s1mat <- as.matrix(fithierncp, pars = "S1") 
z2mat <- as.matrix(fithierncp, pars = "Z2") 
prmat <- as.matrix(fithierncp, pars = "Pr")
logsmoltsmat <- as.matrix(fithierncp, pars = "logsmolts_true")

# estimating mu1 and mu2 but in matrix form
mu1mat <- logsmoltsmat + log(prmat) - z1mat
mu2mat <- logsmoltsmat + log(1 - prmat) - z1mat - z2mat

z1tib <- z1mat %>% t %>% as_tibble %>% mutate(river_name = river_names)
z2tib <- z2mat %>% t %>% as_tibble %>% mutate(river_name = river_names)
prtib <-  prmat %>% t %>% as_tibble %>%  mutate(river_name = river_names)
logsmoltstib <- logsmoltsmat %>% t %>% as_tibble %>% mutate(river_name = river_names)
mu1tib <- mu1mat %>% t %>% as_tibble %>% mutate(river_name = river_names)
mu2tib <- mu2mat %>% t %>% as_tibble %>% mutate(river_name = river_names)

## mapply is WAY faster
# microbenchmark::microbenchmark(diag(cor(z1mat, mu1mat)), 
#                               mapply(cor, as.data.frame(z1mat), as.data.frame(mu1mat)), 
#                               times = 20)

mu1out <- tibble(river_name = rep(unique(river_names), each = 4),
                 correlation = rep(c("mu Z1", "mu Pr", "mu logsmolts", "mu sum"), 7),
       median = NaN,
       lowerCI = NaN,
       upperCI = NaN,
       meansum = NaN)

mu2out <- tibble(river_name = rep(unique(river_names), each = 5),
                 correlation = rep(c("mu2 Z1", "mu2 Z2", "mu2 Pr", "mu2 logsmolts", "mu2 sum"), 7),
                 median = NaN,
                 lowerCI = NaN,
                 upperCI = NaN)


# Estimating correlation squared among parameters 
for (i in unique(river_names)) {
  z1sub <- filter(z1tib, river_name == i) %>% dplyr::select(-river_name) 
  prsub <- filter(prtib, river_name == i) %>% dplyr::select(-river_name) 
  logsmoltssub <- filter(logsmoltstib, river_name == i) %>% dplyr::select(-river_name)  
  mu1sub <- filter(mu1tib, river_name == i) %>% dplyr::select(-river_name) 
 #     print(class(z1sub))
  # print(head(mu1sub))
  mu1z1 <- mapply(cor, mu1sub, z1sub)^2 
  mu1pr <- mapply(cor, mu1sub, log(prsub))^2
  mu1logsmolts <-  mapply(cor, mu1sub,logsmoltssub)^2
  rsqsum <- colSums(rbind(mu1z1, mu1pr, mu1logsmolts))
  print(i)
  print(paste("R squared mu with Z1 =", median(mu1z1)))
  print(paste("R squared mu with log(Pr) =", median(mu1pr)))
  print(paste("R squared mu with logsmolts_true =", median(mu1logsmolts)))
  print(paste("median R squared sum =", median(rsqsum)))
  mu1out[mu1out$river_name == i & mu1out$correlation == "mu Z1", 3] <- median(mu1z1)
  mu1out[mu1out$river_name == i & mu1out$correlation == "mu Z1", 4] <- quantile(mu1z1, 0.05)
  mu1out[mu1out$river_name == i & mu1out$correlation == "mu Z1", 5] <- quantile(mu1z1, 0.95)
  mu1out[mu1out$river_name == i & mu1out$correlation == "mu Pr", 3] <- median(mu1pr)
  mu1out[mu1out$river_name == i & mu1out$correlation == "mu Pr", 4] <- quantile(mu1pr, 0.05)
  mu1out[mu1out$river_name == i & mu1out$correlation == "mu Pr", 5] <- quantile(mu1pr, 0.95)  
  mu1out[mu1out$river_name == i & mu1out$correlation == "mu logsmolts", 3] <- median(mu1logsmolts)
  mu1out[mu1out$river_name == i & mu1out$correlation == "mu logsmolts", 4] <- quantile(mu1logsmolts, 0.05)
  mu1out[mu1out$river_name == i & mu1out$correlation == "mu logsmolts", 5] <- quantile(mu1logsmolts, 0.95)
  mu1out[mu1out$river_name == i & mu1out$correlation == "mu sum", 3] <- median(rsqsum)
  mu1out[mu1out$river_name == i & mu1out$correlation == "mu sum", 4] <- quantile(rsqsum, 0.05)
  mu1out[mu1out$river_name == i & mu1out$correlation == "mu sum", 5] <- quantile(rsqsum, 0.95)  
}


for (i in unique(river_names)) {
  z1sub <- filter(z1tib, river_name == i) %>% select(-river_name)
  z2sub <- filter(z2tib, river_name == i) %>% select(-river_name)
  prsub <- filter(prtib, river_name == i) %>% select(-river_name) 
  logsmoltssub <- filter(logsmoltstib, river_name == i) %>% select(-river_name) 
  mu2sub <- filter(mu2tib, river_name == i) %>% select(-river_name)
 # print(head(mu1sub))
  mu2z1 <- mapply(cor, mu2sub, z1sub)^2 
  mu2z2 <- mapply(cor, mu2sub, z2sub)^2 
  mu2pr <- mapply(cor, mu2sub, log(1 - prsub))^2
  mu2logsmolts <-  mapply(cor, mu2sub, logsmoltssub)^2
  rsqsum2 <- colSums(rbind(mu2z1, mu2z2, mu2pr, mu2logsmolts))
  print(i)
  print(paste("R squared mu with Z1 =", median(mu2z1)))
  print(paste("R squared mu with Z2 =", median(mu2z2)))
  print(paste("R squared mu with log(Pr) =", median(mu2pr)))
  print(paste("R squared mu with logsmolts_true =", median(mu2logsmolts)))
  print(paste("median R squared sum =", median(rsqsum2)))
  mu2out[mu2out$river_name == i & mu2out$correlation == "mu2 Z1", 3] <- median(mu2z1)
  mu2out[mu2out$river_name == i & mu2out$correlation == "mu2 Z1", 4] <- quantile(mu2z1, 0.05)
  mu2out[mu2out$river_name == i & mu2out$correlation == "mu2 Z1", 5] <- quantile(mu2z1, 0.95)  
  mu2out[mu2out$river_name == i & mu2out$correlation == "mu2 Z2", 3] <- median(mu2z2)
  mu2out[mu2out$river_name == i & mu2out$correlation == "mu2 Z2", 4] <- quantile(mu2z2, 0.05)
  mu2out[mu2out$river_name == i & mu2out$correlation == "mu2 Z2", 5] <- quantile(mu2z2, 0.95)
  mu2out[mu2out$river_name == i & mu2out$correlation == "mu2 Pr", 3] <- median(mu2pr)
  mu2out[mu2out$river_name == i & mu2out$correlation == "mu2 Pr", 4] <- quantile(mu2pr, 0.05)
  mu2out[mu2out$river_name == i & mu2out$correlation == "mu2 Pr", 5] <- quantile(mu2pr, 0.95)  
  mu2out[mu2out$river_name == i & mu2out$correlation == "mu2 logsmolts", 3] <- median(mu2logsmolts)
  mu2out[mu2out$river_name == i & mu2out$correlation == "mu2 logsmolts", 4] <- quantile(mu2logsmolts, 0.05)
  mu2out[mu2out$river_name == i & mu2out$correlation == "mu2 logsmolts", 5] <- quantile(mu2logsmolts, 0.95)
  mu2out[mu2out$river_name == i & mu2out$correlation == "mu2 sum", 3] <- median(rsqsum2)
  mu2out[mu2out$river_name == i & mu2out$correlation == "mu2 sum", 4] <- quantile(rsqsum2, 0.05)
  mu2out[mu2out$river_name == i & mu2out$correlation == "mu2 sum", 5] <- quantile(rsqsum2, 0.95)  
}

### PLOTTING PSEUDO R^2s

print(mu1out, n = 28)
print(mu2out, n = 28)

mu1out %>%
  mutate(param = word(correlation, 2)) %>%
  mutate(river_name = fct_relevel(river_name, propgrilse$river_name)) %>%
  filter(correlation != "mu sum") %>%
  ggplot(aes(param, median, color = river_name)) + 
  geom_point(position = position_dodge(-0.4)) + 
  geom_errorbar(aes(ymax = upperCI, ymin = lowerCI), width = 0.0, position = position_dodge(-0.4)) +
  ylab("R squared estimates") + xlab("Parameter") + ggtitle("Correlation with R_1") +
  coord_flip() + theme_cowplot()


mu2out %>%
  mutate(param = word(correlation, 2)) %>%
  mutate(river_name = fct_relevel(river_name, propgrilse$river_name)) %>%
  filter(correlation != "mu sum") %>%
  filter(!river_name %in% c("Campbellton River", "Conne River", "Western Arm Brook")) %>%
  ggplot(aes(param, median, color = river_name)) + 
  geom_point(position = position_dodge(-0.4)) + 
  geom_errorbar(aes(ymax = upperCI, ymin = lowerCI), width = 0.0, position = position_dodge(-0.4)) +
  ylab("R squared estimates") + xlab("Parameter") + ggtitle("Correlation with R_2") +
  coord_flip() + theme_cowplot()

param_names <- c("mu" = "italic(R)['1,est']", "mu2" = "italic(R)['2,est']")

bind_rows(mu1out, mu2out) %>%
  mutate(ret = word(correlation, 1), param = word(correlation, 2)) %>%
  mutate(river_name = fct_relevel(river_name, propgrilse$river_name)) %>%
  mutate(param = forcats::fct_relevel(param, rev(c("logsmolts", "Z1", "Z2", "Pr", "sum")))) %>%
  # Need to replace with NA rather than filter out so the rivers align in both facets
  mutate(median = ifelse((river_name %in% c("Campbellton River", "Conne River", "Western Arm Brook") & ret == "mu2"),
NA, median)) %>%
  mutate(lowerCI = ifelse((river_name %in% c("Campbellton River", "Conne River", "Western Arm Brook") & ret == "mu2"),
                         NA, lowerCI)) %>%
  mutate(upperCI = ifelse((river_name %in% c("Campbellton River", "Conne River", "Western Arm Brook") & ret == "mu2"),
                         NA, upperCI)) %>%
  #filter(!(river_name %in% c("Campbellton River", "Conne River", "Western Arm Brook") & ret == "mu2")) %>%
  filter(correlation != "mu sum" ) %>%
  filter(correlation != "mu2 sum" ) %>%
  ggplot(aes(param, median, color = river_name)) + 
  geom_point(size = 2.5, position = position_dodge(-0.9)) +
  #geom_vline(aes(xintercept = 1.5), color = "gray70", alpha = 0.5) +
  geom_vline(aes(xintercept = 1.5), color = "gray70", alpha = 0.5) +
  geom_vline(aes(xintercept = 2.5), color = "gray70", alpha = 0.5) +
  geom_vline(aes(xintercept = 3.5), color = "gray70", alpha = 0.5) +
#  geom_vline(aes(xintercept = 4.5), color = "gray70", alpha = 0.5) +
  geom_hline(aes(yintercept = 1), color = "gray70", alpha = 0.5) +
  geom_hline(aes(yintercept = 0), color = "gray70", alpha = 0.5) +
  geom_errorbar(aes(ymax = upperCI, ymin = lowerCI), width = 0.0, size = 0.7, position = position_dodge(-0.9)) +
  ylab(bquote("Estimated R"^2)) + 
  xlab("Parameter") + labs(color = "Population") +
  scale_x_discrete(labels = c('logsmolts' = expression(italic(log(smolts))),
                              'Z1'= expression(italic(Z)[1]),
                              'Z2'= expression(italic(Z)[2]),
                              'Pr'= expression(italic(P)[g]))) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                     breaks = seq(0,1, by = 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_flip() + facet_grid(~ret, scales = "free_x", 
                            labeller = as_labeller(param_names, label_parsed)) + theme_cowplot()
ggsave("figures/rsquared-mu-params.png", width = 7.5, height = 5) 







    #########################


mu1var <- tibble(river_name = rep(unique(river_names), each = 7),
                 correlation = rep(c("Var(logsmolts)/Var(mu)", 
                                     "Var(Z1)/Var(mu)", 
                                     "Var(logPr)/Var(mu)",
                                     "2*Cov(logsmolts,-Z1)/Var(mu)",
                                     "2*Cov(logsmolts,logPr)/Var(mu)",
                                     "2*Cov(-Z1,logPr)/Var(mu)",
                                     "mu sum"), 7),
                 median = NaN,
                 lowerCI = NaN,
                 upperCI = NaN)

# calculating sum of correlated variables
for (i in unique(river_names)) {
  z1sub <- filter(z1tib, river_name == i) %>% dplyr::select(-river_name) 
  prsub <- filter(prtib, river_name == i) %>% dplyr::select(-river_name) 
  logsmoltssub <- filter(logsmoltstib, river_name == i) %>% dplyr::select(-river_name)  
  mu1sub <- filter(mu1tib, river_name == i) %>% dplyr::select(-river_name) 
  #     print(class(z1sub))
  # print(head(mu1sub))
  varmu1 <- mapply(var, mu1sub)
  varlogpr <- mapply(var, log(prsub))
  varz1 <- mapply(var, z1sub)
  varlogsmolts <- mapply(var, logsmoltssub)
  covlogsmoltsz1 <- mapply(cov, logsmoltssub, z1sub*(-1)) 
  covlogsmoltslogpr <- mapply(cov, logsmoltssub, log(prsub)) 
  covz1logpr <- mapply(cov, z1sub*(-1),  log(prsub)) 
  
  varsmoltsmu <- varlogsmolts/varmu1
  varz1mu <- varz1/varmu1
  varlogprmu <- varlogpr/varmu1
  cov2logsmoltsz1mu <- 2*covlogsmoltsz1/varmu1
  cov2logsmoltslogprmu <- 2*covlogsmoltslogpr/varmu1
  cov2z1logprmu <- 2*covz1logpr/varmu1
  
  rsqsum <- colSums(rbind(varsmoltsmu, varz1mu, varlogprmu, cov2logsmoltsz1mu, cov2logsmoltslogprmu, cov2z1logprmu))
  
  print(i)
  # print(paste("R squared mu with Z1 =", median(mu1z1)))
  # print(paste("R squared mu with log(Pr) =", median(mu1pr)))
  # print(paste("R squared mu with logsmolts_true =", median(mu1logsmolts)))
  # print(paste("median R squared sum =", median(rsqsum)))
  mu1var[mu1var$river_name == i & mu1var$correlation == "Var(logsmolts)/Var(mu)", 3] <- median(varsmoltsmu)
  mu1var[mu1var$river_name == i & mu1var$correlation == "Var(logsmolts)/Var(mu)", 4] <- quantile(varsmoltsmu, 0.05)
  mu1var[mu1var$river_name == i & mu1var$correlation == "Var(logsmolts)/Var(mu)", 5] <- quantile(varsmoltsmu, 0.95)
  mu1var[mu1var$river_name == i & mu1var$correlation == "Var(Z1)/Var(mu)", 3] <- median(varz1mu)
  mu1var[mu1var$river_name == i & mu1var$correlation == "Var(Z1)/Var(mu)", 4] <- quantile(varz1mu, 0.05)
  mu1var[mu1var$river_name == i & mu1var$correlation == "Var(Z1)/Var(mu)", 5] <- quantile(varz1mu, 0.95)  
  mu1var[mu1var$river_name == i & mu1var$correlation == "Var(logPr)/Var(mu)", 3] <- median(varlogprmu)
  mu1var[mu1var$river_name == i & mu1var$correlation == "Var(logPr)/Var(mu)", 4] <- quantile(varlogprmu, 0.05)
  mu1var[mu1var$river_name == i & mu1var$correlation == "Var(logPr)/Var(mu)", 5] <- quantile(varlogprmu, 0.95)
  mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(logsmolts,-Z1)/Var(mu)", 3] <- median(cov2logsmoltsz1mu)
  mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(logsmolts,-Z1)/Var(mu)", 4] <- quantile(cov2logsmoltsz1mu, 0.05)
  mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(logsmolts,-Z1)/Var(mu)", 5] <- quantile(cov2logsmoltsz1mu, 0.95)
  mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(logsmolts,logPr)/Var(mu)", 3] <- median(cov2logsmoltslogprmu)
  mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(logsmolts,logPr)/Var(mu)", 4] <- quantile(cov2logsmoltslogprmu, 0.05)
  mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(logsmolts,logPr)/Var(mu)", 5] <- quantile(cov2logsmoltslogprmu, 0.95)
  mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(-Z1,logPr)/Var(mu)", 3] <- median(cov2z1logprmu)
  mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(-Z1,logPr)/Var(mu)", 4] <- quantile(cov2z1logprmu, 0.05)
  mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(-Z1,logPr)/Var(mu)", 5] <- quantile(cov2z1logprmu, 0.95)
  mu1var[mu1var$river_name == i & mu1var$correlation == "mu sum", 3] <- median(rsqsum)
  mu1var[mu1var$river_name == i & mu1var$correlation == "mu sum", 4] <- quantile(rsqsum, 0.05)
  mu1var[mu1var$river_name == i & mu1var$correlation == "mu sum", 5] <- quantile(rsqsum, 0.95)  
}

mu1var

mu1var %>%
  mutate(river_name = fct_relevel(river_name, propgrilse$river_name),
         correlation = fct_relevel(correlation, "mu sum", after = Inf)) %>%
  ggplot(aes(correlation, median, color = river_name)) + 
  geom_point(position = position_dodge(-0.5)) + 
  geom_errorbar(aes(ymax = upperCI, ymin = lowerCI), width = 0.0, position = position_dodge(-0.5)) +
  ylab("Contribution to Variance") + xlab("Parameter") + ggtitle("Correlation with mu_1") +
  coord_flip() + theme_cowplot()
ggsave("figures/correlatedvariancemu1.png", height = 4.5, width = 8)


### 
# mu1mat <- logsmoltsmat + log(prmat) - z1mat
# mu2mat <- logsmoltsmat + log(1 - prmat) - z1mat - z2mat

#Var(mu) = Var(logsmolts) + Var(Z1) + Var(logPr) + 2*Cov(logsmolts, -Z1) + 2*Cov(logsmolts, logPr) + 2*Cov(-Z1, logPr)
#Var(mu2) = Var(logsmolts) + Var(Z1) + Var(Z2) + Var(log(1-Pr)) + 
#           2*Cov(logsmolts, -Z1) + 2*Cov(logsmolts, log(1-Pr)) + 2*Cov(logsmolts, -Z2) +
#           2*Cov(-Z1,  -Z2) + 2*Cov(-Z1,  log(1-Pr)) +   2*Cov(-Z2,  log(1-Pr))



mu2var <- tibble(river_name = rep(unique(river_names), each = 11),
                 correlation = rep(c("Var(logsmolts)/Var(mu2)", 
                                     "Var(Z1)/Var(mu2)", 
                                     "Var(Z2)/Var(mu2)", 
                                     "Var(log(1-Pr))/Var(mu2)",
                                     "2*Cov(logsmolts,-Z1)/Var(mu2)",
                                     "2*Cov(logsmolts,-Z2)/Var(mu2)",
                                     "2*Cov(logsmolts,log(1-Pr))/Var(mu2)",
                                     "2*Cov(-Z1,log(1-Pr))/Var(mu2)",
                                     "2*Cov(-Z2,log(1-Pr))/Var(mu2)",
                                     "2*Cov(-Z1,-Z2)/Var(mu2)",
                                     "mu2 sum"), 7),
                 median = NaN,
                 lowerCI = NaN,
                 upperCI = NaN)

# calculating sum of correlated variables
for (i in unique(river_names)) {
  z1sub <- filter(z1tib, river_name == i) %>% dplyr::select(-river_name) 
  z2sub <- filter(z2tib, river_name == i) %>% dplyr::select(-river_name) 
  prsub <- filter(prtib, river_name == i) %>% dplyr::select(-river_name) 
  logsmoltssub <- filter(logsmoltstib, river_name == i) %>% dplyr::select(-river_name)  
  mu2sub <- filter(mu2tib, river_name == i) %>% dplyr::select(-river_name) 
  #     print(class(z1sub))
  # print(head(mu1sub))
  varmu2 <- mapply(var, mu2sub)
  
  varlogpr <- mapply(var, log(1-prsub))
  varz1 <- mapply(var, z1sub)
  varz2 <- mapply(var, z2sub)
  varlogsmolts <- mapply(var, logsmoltssub)
  
  covlogsmoltsz1 <- mapply(cov, logsmoltssub, z1sub*(-1)) 
  covlogsmoltsz2 <- mapply(cov, logsmoltssub, z2sub*(-1)) 
  covlogsmoltslogpr <- mapply(cov, logsmoltssub, log(1-prsub)) 
  covz1logpr <- mapply(cov, z1sub*(-1),  log(1-prsub)) 
  covz2logpr <- mapply(cov, z2sub*(-1),  log(1-prsub)) 
  covz1z2 <- mapply(cov, z1sub*(-1),  z2sub*(-1)) 
  
  varsmoltsmu <- varlogsmolts/varmu2
  varz1mu <- varz1/varmu2
  varz2mu <- varz2/varmu2
  varlogprmu <- varlogpr/varmu2
  
  cov2logsmoltsz1mu <- 2*covlogsmoltsz1/varmu2
  cov2logsmoltsz2mu <- 2*covlogsmoltsz2/varmu2
  cov2logsmoltslogprmu <- 2*covlogsmoltslogpr/varmu2
  cov2z1logprmu <- 2*covz1logpr/varmu2
  cov2z2logprmu <- 2*covz2logpr/varmu2
  cov2z1z2 <- 2*covz1z2/varmu2
  
  rsqsum <- colSums(rbind(varsmoltsmu, varz1mu, varz2mu, varlogprmu, 
                          cov2logsmoltsz1mu, cov2logsmoltsz2mu, cov2logsmoltslogprmu, 
                          cov2z1logprmu, cov2z2logprmu, cov2z1z2))
  #range(rsqsum)
  print(i)
  # print(paste("R squared mu with Z1 =", median(mu1z1)))
  # print(paste("R squared mu with log(Pr) =", median(mu1pr)))
  # print(paste("R squared mu with logsmolts_true =", median(mu1logsmolts)))
  # print(paste("median R squared sum =", median(rsqsum)))
  mu2var[mu2var$river_name == i & mu2var$correlation == "Var(logsmolts)/Var(mu2)", 3] <- median(varsmoltsmu)
  mu2var[mu2var$river_name == i & mu2var$correlation == "Var(logsmolts)/Var(mu2)", 4] <- quantile(varsmoltsmu, 0.05)
  mu2var[mu2var$river_name == i & mu2var$correlation == "Var(logsmolts)/Var(mu2)", 5] <- quantile(varsmoltsmu, 0.95)
  mu2var[mu2var$river_name == i & mu2var$correlation == "Var(Z1)/Var(mu2)", 3] <- median(varz1mu)
  mu2var[mu2var$river_name == i & mu2var$correlation == "Var(Z1)/Var(mu2)", 4] <- quantile(varz1mu, 0.05)
  mu2var[mu2var$river_name == i & mu2var$correlation == "Var(Z1)/Var(mu2)", 5] <- quantile(varz1mu, 0.95)  
  mu2var[mu2var$river_name == i & mu2var$correlation == "Var(Z2)/Var(mu2)", 3] <- median(varz2mu)
  mu2var[mu2var$river_name == i & mu2var$correlation == "Var(Z2)/Var(mu2)", 4] <- quantile(varz2mu, 0.05)
  mu2var[mu2var$river_name == i & mu2var$correlation == "Var(Z2)/Var(mu2)", 5] <- quantile(varz2mu, 0.95)  
  mu2var[mu2var$river_name == i & mu2var$correlation == "Var(log(1-Pr))/Var(mu2)", 3] <- median(varlogprmu)
  mu2var[mu2var$river_name == i & mu2var$correlation == "Var(log(1-Pr))/Var(mu2)", 4] <- quantile(varlogprmu, 0.05)
  mu2var[mu2var$river_name == i & mu2var$correlation == "Var(log(1-Pr))/Var(mu2)", 5] <- quantile(varlogprmu, 0.95)
  
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(logsmolts,-Z1)/Var(mu2)", 3] <- median(cov2logsmoltsz1mu)
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(logsmolts,-Z1)/Var(mu2)", 4] <- quantile(cov2logsmoltsz1mu, 0.05)
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(logsmolts,-Z1)/Var(mu2)", 5] <- quantile(cov2logsmoltsz1mu, 0.95)
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(logsmolts,-Z2)/Var(mu2)", 3] <- median(cov2logsmoltsz2mu)
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(logsmolts,-Z2)/Var(mu2)", 4] <- quantile(cov2logsmoltsz2mu, 0.05)
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(logsmolts,-Z2)/Var(mu2)", 5] <- quantile(cov2logsmoltsz2mu, 0.95)  
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(logsmolts,log(1-Pr))/Var(mu2)", 3] <- median(cov2logsmoltslogprmu)
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(logsmolts,log(1-Pr))/Var(mu2)", 4] <- quantile(cov2logsmoltslogprmu, 0.05)
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(logsmolts,log(1-Pr))/Var(mu2)", 5] <- quantile(cov2logsmoltslogprmu, 0.95)
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(-Z1,log(1-Pr))/Var(mu2)", 3] <- median(cov2z1logprmu)
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(-Z1,log(1-Pr))/Var(mu2)", 4] <- quantile(cov2z1logprmu, 0.05)
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(-Z1,log(1-Pr))/Var(mu2)", 5] <- quantile(cov2z1logprmu, 0.95)  
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(-Z2,log(1-Pr))/Var(mu2)", 3] <- median(cov2z2logprmu)
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(-Z2,log(1-Pr))/Var(mu2)", 4] <- quantile(cov2z2logprmu, 0.05)
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(-Z2,log(1-Pr))/Var(mu2)", 5] <- quantile(cov2z2logprmu, 0.95)
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(-Z1,-Z2)/Var(mu2)", 3] <- median(cov2z1z2)
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(-Z1,-Z2)/Var(mu2)", 4] <- quantile(cov2z1z2, 0.05)
  mu2var[mu2var$river_name == i & mu2var$correlation == "2*Cov(-Z1,-Z2)/Var(mu2)", 5] <- quantile(cov2z1z2, 0.95)
  
  mu2var[mu2var$river_name == i & mu2var$correlation == "mu2 sum", 3] <- median(rsqsum)
  mu2var[mu2var$river_name == i & mu2var$correlation == "mu2 sum", 4] <- quantile(rsqsum, 0.05)
  mu2var[mu2var$river_name == i & mu2var$correlation == "mu2 sum", 5] <- quantile(rsqsum, 0.95)  
}

print(mu2var, n = 11)

mu2var %>%
  mutate(river_name = fct_relevel(river_name, propgrilse$river_name),
         correlation = fct_relevel(correlation, "mu2 sum", after = Inf)) %>%
  ggplot(aes(correlation, median, color = river_name)) + 
  geom_point(position = position_dodge(-0.5)) + 
  geom_errorbar(aes(ymax = upperCI, ymin = lowerCI), width = 0.0, position = position_dodge(-0.5)) +
  ylab("Contribution to Variance") + xlab("Parameter") + ggtitle("Correlation with mu_2") +
  coord_flip() + theme_cowplot()
ggsave("figures/correlatedvariancemu2.png", height = 6, width = 8)

mu1var


saveRDS(mu1var, file = "data/mu1-variance-decomposed.rds")
saveRDS(mu2var, file = "data/mu2-variance-decomposed.rds")




rsqspost <- z1post %>%
  rename(z1 = value) %>%
  mutate(z2 = z2post$value,    # joining posteriors
         pr = prpost$value, 
         logsmolts = logsmoltspost$value,
         mu1 = mu1, mu2 = mu2) %>%
  group_by(river_name, year) %>%
  summarise(cor1z1 = cor(mu1, z1)^2,   # estimating R^2s
            cor1pr = cor(mu1, log(pr))^2,
            cor1smolts = cor(mu1, logsmolts)^2,
            cor2z1 = cor(mu2, z1)^2,
            cor2z2 = cor(mu2, z2)^2,
            cor2pr = cor(mu2, log(1-pr))^2,
            cor2smolts = cor(mu2, logsmolts)^2)

rsqspost %>%
  ungroup %>% group_by(river_name) %>%
  summarise_at(vars(starts_with("cor")), mean) # calculating mean R^2 among years

rsqspost %>%
  ungroup %>% group_by(river_name) %>%
  summarise_at(vars(starts_with("cor")), quantile, c(0.05)) # calculating upper 90% CI

rsqspost %>%
ungroup %>% group_by(river_name) %>%
  summarise_at(vars(starts_with("cor")), quantile, c(0.95)) # calculating lower 90% CI

# As above but with empirical log return estimates from scale data (rather than mu and mu2 from model equations)
rsqspost2 <- z1post %>%
  rename(z1 = value) %>%
  mutate(z2 = z2post$value,    # joining posteriors
         pr = prpost$value, 
         logsmolts = logsmoltspost$value) %>%
  left_join(estreturns, by = c("year", "river_name"))

rsqspost2 %>%
  group_by(river_name) %>%
  summarise(min = min(log2SW), max = max(log2SW))

rsqspost2 %>%
  group_by(river_name, year) %>%
  summarise(cor1z1 = cor(log1SW, z1)^2,   # estimating R^2s
            cor1pr = cor(log1SW, log(pr))^2,
            cor1smolts = cor(log1SW, logsmolts)^2,
            cor2z1 = cor(log2SW, z1)^2,
            cor2z2 = cor(log2SW, z2)^2,
            cor2pr = cor(log2SW, log(1-pr))^2,
            cor2smolts = cor(log2SW, logsmolts)^2)


rsqspost2 %>%
  ungroup %>% group_by(river_name) %>%
  summarise_at(vars(starts_with("cor")), mean) # calculating mean R^2 among years

rsqspost2 %>%
  ungroup %>% group_by(river_name) %>%
  summarise_at(vars(starts_with("cor")), quantile, c(0.05)) # calculating upper 90% CI

rsqspost2 %>%
  ungroup %>% group_by(river_name) %>%
  summarise_at(vars(starts_with("cor")), quantile, c(0.95)) # calculating lower 90% CI

######


river_names <- s1quant %>%
  ungroup() %>%
  select(year, zscore, river_name) %>%
  spread(year, zscore)  %>%
  pull(river_name)

s1mat <- s1quant %>%
  #filter(year <= 2013) %>%
  ungroup() %>%
  select(year, zscore, river_name) %>%
  spread(year, zscore) %>%
  select(-river_name) %>%
  t %>% as.matrix

colnames(s1mat) <- c(word(river_names)[1:6], "WAB")

M <- cor(s1mat, use = "complete.obs", method = c("pearson"))
M

# matrix of the p-value of the correlation
# Pearson's product moment correlation coefficient
p.mat <- cor.mtest(s1mat)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

oldpar <- par() 

png("figures/corr-s1.png", width = 7.5, height = 3, units = "in", res = 300, pointsize = 9.5)
par(mfrow = c(1,2), mar = c(0,0,0,0), oma = c(0,0,0,0))
corrplot(M, method = "color", col = col(200),  
         type = "upper", order = "hclust",  diag = FALSE, 
          number.cex = 0.8,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", # only coloured boxes have significant pairwise p-values.
         # hide correlation coefficient on the principal diagonal
         # diag=FALSE,
       #  title = expression(paste("Correlation of standardized Z-scores in ", S[1])),
         mar = c(0, 0, 0, 0)
)

mtext("(a)", line = -2, adj =0.1, cex = 1.6)

corrplot(M, method = "color", col = col(200),  
         type = "upper", order = "hclust",  diag = FALSE, 
         number.cex = 0.8,
         #addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "label_sig", # only coloured boxes have significant pairwise p-values.
         # hide correlation coefficient on the principal diagonal
         # diag=FALSE,
         #title = expression(paste("Correlation of standardized Z-scores in ", S[1])),
         mar = c(0, 0, 0, 0)
)
  mtext("(b)", line = -2, adj =0.1, cex = 1.6)
dev.off()
    

ggplot(s1quant, aes(year, median, color = river_name)) + 
  geom_line(alpha = 0.5) + geom_point(size = 2) + ylab("S1 estimates") +
  geom_errorbar(aes(ymin = q05, ymax = q95),width= 0.1)


ggplot(s1quant, aes(year, median, color = river_name)) + 
  geom_line() + geom_point(size = 2) + ylab("S1 estimates") +
  geom_errorbar(aes(ymin = q05, ymax = q95),width= 0.1) + 
  facet_wrap(~river_name)


### Correlation of smolts and returns

smoltsmat <- allreturns %>%
  select(river_name, year, zscoresmolts) %>%
  spread(year, zscoresmolts) %>%
  ungroup() %>%
  select(-river_name) %>%
  t %>% as.matrix

totalmat <- allreturns %>%
  select(river_name, year, zscoretotal) %>%
  spread(year, zscoretotal) %>%
  ungroup() %>%
  select(-river_name) %>%
  t %>% as.matrix



colnames(smoltsmat) <- river_names
colnames(totalmat) <- river_names


Msmolts <- cor(smoltsmat, use = "complete.obs", method = c("pearson"))
Msmolts
Mtotal <- cor(totalmat, use = "complete.obs", method = c("pearson"))
Mtotal





# matrix of the p-value of the correlation
p.mat.smolts <- cor.mtest(smoltsmat)
p.mat.total <- cor.mtest(totalmat)

tm <- cor.test(smoltsmat[,4], smoltsmat[,5])


col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(Msmolts, method = "color", col = col(200),  
         type = "upper", order = "hclust",  diag = FALSE,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         # p.mat = p.mat.smolts, sig.level = 0.01, insig = "blank", # only coloured boxes have significant pairwise p-values.
         # hide correlation coefficient on the principal diagonal
         # diag=FALSE,
         mar = c(0, 0, 1.8, 0),
         title = "Correlation of standardized Z-scores of log-transformed smolt abundances"
)

corrplot(Mtotal, method = "color", col = col(200),  
         type = "upper", order = "hclust", diag = FALSE,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat.total, sig.level = 0.01, insig = "blank", # only coloured boxes have significant pairwise p-values.
         # hide correlation coefficient on the principal diagonal
         # diag=FALSE,
         mar = c(0, 0, 1.8, 0),
         title = "Correlation of standardized Z-scores of log-transformed total return abundances"
)

allreturns  %>%
  group_by(river_name) %>%
  mutate(zscoresmolts = scale(smolts), zscoretotal = scale(totalreturns)) %>%
  ggplot(aes(year, smolts, colour = river_name)) + geom_point() + geom_line() +
  facet_wrap(~river_name, scale = "free_y")

allreturns  %>%
  group_by(river_name) %>%
  mutate(zscoresmolts = scale(smolts), zscoretotal = scale(totalreturns)) %>%
  ggplot(aes(year, zscoresmolts, colour = river_name)) + geom_point() + geom_line()

