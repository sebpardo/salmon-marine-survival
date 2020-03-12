
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

mu1errorvar <- tibble(river_name = rep(unique(river_names), each = 7),
                      correlation = rep(c("Var(logsmolts)", 
                                          "Var(Z1)", 
                                          "Var(logPr)",
                                          "Var(mu)",
                                          "Cov(logsmolts,-Z1)",
                                          "Cov(logsmolts,logPr)",
                                          "Cov(-Z1,logPr)"), 7),
                      mean = NaN)

mu1correctedcov <- tibble(river_name = rep(unique(river_names), each = 3),
                 correlation = rep(c("Corrected cov(logsmolts,-Z1)",
                                     "Corrected cov(logsmolts,logPr)",
                                     "Corrected cov(-Z1,logPr)"), 7),
                 median = NaN,
                 lowerCI = NaN,
                 upperCI = NaN)

mu1correctedrsq <- tibble(river_name = rep(unique(river_names), each = 4),
                          correlation = rep(c("R2(mu1,-Z1)",
                                              "R2(mu1,logsmolts)",
                                              "R2(mu1,logPr)",
                                              "R2 sum"), 7),
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
  
  covmu1logsmolts <- mapply(cov, mu1sub, logsmoltssub) 
  covmu1z1 <- mapply(cov, mu1sub, z1sub*(-1)) 
  covmu1logpr <- mapply(cov, mu1sub,  log(prsub)) 
  
  
  # estimating error variances and covariances
  # this is done by estimating variance across iterations for each year and then averaging across years,
  # hence the transpose
  errorvarmu1 <- mapply(var, mu1sub %>% t %>% as_tibble) %>% mean
  errorvarlogpr <- mapply(var, log(prsub) %>% t %>% as_tibble) %>% mean
  errorvarz1 <- mapply(var, z1sub %>% t %>% as_tibble) %>% mean
  errorvarlogsmolts <- mapply(var, logsmoltssub %>% t %>% as_tibble) %>% mean
  
  errorcovlogsmoltsz1 <- mapply(cov, logsmoltssub %>% t %>% as_tibble, t(z1sub*(-1)) %>% as_tibble)  %>% mean
  errorcovlogsmoltslogpr <- mapply(cov, logsmoltssub %>% t %>% as_tibble, log(prsub) %>% t %>% as_tibble)  %>% mean
  errorcovz1logpr <- mapply(cov, t(z1sub*(-1)) %>% as_tibble,  log(prsub) %>% t %>% as_tibble)  %>% mean
  
  errorcovmu1logsmolts <- mapply(cov, mu1sub %>% t %>% as_tibble,logsmoltssub %>% t %>% as_tibble)  %>% mean
  errorcovmu1z1 <- mapply(cov, mu1sub %>% t %>% as_tibble, t(z1sub*(-1)) %>% as_tibble)  %>% mean
  errorcovmu1logpr <- mapply(cov, mu1sub %>% t %>% as_tibble, log(prsub) %>% t %>% as_tibble)  %>% mean
  
  # saving to out
  mu1errorvar[mu1errorvar$river_name == i & mu1errorvar$correlation == "Var(mu)", 3] <- errorvarmu1
  mu1errorvar[mu1errorvar$river_name == i & mu1errorvar$correlation == "Var(logPr)", 3] <- errorvarlogpr
  mu1errorvar[mu1errorvar$river_name == i & mu1errorvar$correlation == "Var(Z1)", 3] <- errorvarz1
  mu1errorvar[mu1errorvar$river_name == i & mu1errorvar$correlation == "Var(logsmolts)", 3] <- errorvarlogsmolts
  
  mu1errorvar[mu1errorvar$river_name == i & mu1errorvar$correlation == "Cor(logsmolts,-Z1)", 3] <- errorcovlogsmoltsz1
  mu1errorvar[mu1errorvar$river_name == i & mu1errorvar$correlation == "Cor(logsmolts,logPr)", 3] <- errorcovlogsmoltslogpr
  mu1errorvar[mu1errorvar$river_name == i & mu1errorvar$correlation == "Cor(-Z1,logPr)", 3] <- errorcovz1logpr
  
  mu1errorvar[mu1errorvar$river_name == i & mu1errorvar$correlation == "Cor(mu1, logsmolts)", 3] <- errorcovmu1logsmolts
  mu1errorvar[mu1errorvar$river_name == i & mu1errorvar$correlation == "Cor(mu1, -Z1)", 3] <- errorcovmu1z1
  mu1errorvar[mu1errorvar$river_name == i & mu1errorvar$correlation == "Cor(mu1, logPr)", 3] <- errorcovmu1logpr
  
  # estimating corrected correlation values
  corrected_corlogsmoltsz1 <- (covlogsmoltsz1 - errorcovlogsmoltsz1)/(sqrt(varlogsmolts - errorvarlogsmolts)*sqrt(varz1 - errorvarz1))
  corrected_corlogsmoltslogpr <- (covlogsmoltslogpr - errorcovlogsmoltslogpr)/(sqrt(varlogsmolts - errorvarlogsmolts)*sqrt(varlogpr - errorvarlogpr))
  corrected_corz1logpr <- (covz1logpr - errorcovz1logpr)/(sqrt(varlogpr - errorvarlogpr)*sqrt(varz1 - errorvarz1))
  
  # estimating R^2 with mu1
  corrected_r2_mu1logsmolts <- ((covmu1logsmolts - errorcovmu1logsmolts)/(sqrt(varlogsmolts - errorvarlogsmolts)*sqrt(varmu1 - errorvarmu1)))^2
  corrected_r2_mu1z1 <- ((covmu1z1 - errorcovmu1z1)/(sqrt(varz1 - errorvarz1)*sqrt(varmu1 - errorvarmu1)))^2
  corrected_r2_mu1logpr <- ((covmu1logpr - errorcovmu1logpr)/(sqrt(varlogpr - errorvarlogpr)*sqrt(varmu1 - errorvarmu1)))^2
  
  rsqsum <- colSums(rbind(corrected_r2_mu1logsmolts, corrected_r2_mu1z1, corrected_r2_mu1logpr))
  
  
  # storing values, na.rm = TRUE to ignore NaNs created when square-rooting a few negative values
  mu1correctedcov[mu1correctedcov$river_name == i & 
                    mu1correctedcov$correlation == "Corrected cov(logsmolts,-Z1)", 3:5] <- quantile(corrected_corlogsmoltsz1, 
                                                                                                    c(0.5, 0.05, 0.95), na.rm = TRUE)
  mu1correctedcov[mu1correctedcov$river_name == i & 
                    mu1correctedcov$correlation == "Corrected cov(logsmolts,logPr)", 3:5] <- quantile(corrected_corlogsmoltslogpr, 
                                                                                                      c(0.5, 0.05, 0.95), na.rm = TRUE)
  mu1correctedcov[mu1correctedcov$river_name == i & 
                    mu1correctedcov$correlation == "Corrected cov(-Z1,logPr)", 3:5] <- quantile(corrected_corz1logpr, 
                                                                                                c(0.5, 0.05, 0.95), na.rm = TRUE)
  
   
    # storing values, na.rm = TRUE to ignore NaNs created when square-rooting a few negative values
  mu1correctedrsq[mu1correctedrsq$river_name == i & 
                    mu1correctedrsq$correlation == "R2(mu1,logsmolts)", 3:5] <- quantile(corrected_r2_mu1logsmolts,
                                                                                         c(0.5, 0.05, 0.95), na.rm = TRUE)
  mu1correctedrsq[mu1correctedrsq$river_name == i & 
                  mu1correctedrsq$correlation == "R2(mu1,-Z1)", 3:5] <- quantile(corrected_r2_mu1z1,
                                                                                 c(0.5, 0.05, 0.95), na.rm = TRUE)
  mu1correctedrsq[mu1correctedrsq$river_name == i & 
                  mu1correctedrsq$correlation == "R2(mu1,logPr)", 3:5] <- quantile(corrected_r2_mu1logpr,
                                                                                   c(0.5, 0.05, 0.95), na.rm = TRUE)
  mu1correctedrsq[mu1correctedrsq$river_name == i & 
                    mu1correctedrsq$correlation == "R2 sum", 3:5] <- quantile(rsqsum, c(0.5, 0.05, 0.95), na.rm = TRUE)
  
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

mu1correctedcov

pivot_wider(mu1errorvar, names_from = correlation, values_from = mean)

mu1correctedcov %>%
  select(-lowerCI, -upperCI) %>%
  mutate(correlation = str_replace(correlation, "cov", "cor")) %>%
  pivot_wider(names_from = correlation, values_from = median)

mu1correctedrsq %>%
  #mutate(correlation = str_replace(correlation, "cov", "cor")) %>%
  mutate(river_name = fct_relevel(river_name, propgrilse$river_name)) %>%
  ggplot(aes(correlation, median, color = river_name)) + 
  geom_hline(yintercept = 1, col = "gray80") +
  geom_point(position = position_dodge(-0.5)) + 
  geom_errorbar(aes(ymax = upperCI, ymin = lowerCI), width = 0.0, position = position_dodge(-0.5)) +
  ylab("R squared value") + xlab("Parameters") + ggtitle("Corrected R squares") +
  labs(color = "River name") +
  coord_flip() + theme_cowplot()
ggsave("figures/corrected_rsq_mu1.png", height = 4.5, width = 8)


mu1correctedcov %>%
  mutate(correlation = str_replace(correlation, "cov", "cor")) %>%
  mutate(river_name = fct_relevel(river_name, propgrilse$river_name)) %>%
  ggplot(aes(correlation, median, color = river_name)) + 
  geom_hline(yintercept = 0, col = "gray80") +
  geom_point(position = position_dodge(-0.5)) + 
  geom_errorbar(aes(ymax = upperCI, ymin = lowerCI), width = 0.0, position = position_dodge(-0.5)) +
  ylab("Correlation value") + xlab("Parameters") + ggtitle("Corrected correlations") +
  labs(color = "River name") + ylim(-1, 1) +
   coord_flip() + theme_cowplot()
ggsave("figures/corrected_cor_mu1.png", height = 4.5, width = 8)


mu1var %>%
  mutate(river_name = fct_relevel(river_name, propgrilse$river_name),
         correlation = fct_relevel(correlation, "mu sum", after = Inf)) %>%
  ggplot(aes(correlation, median, color = river_name)) + 
  geom_hline(yintercept = 0, col = "gray80") +
  geom_point(position = position_dodge(-0.5)) + 
  geom_errorbar(aes(ymax = upperCI, ymin = lowerCI), width = 0.0, position = position_dodge(-0.5)) +
  ylab("Contribution to Variance") + xlab("Parameter") + ggtitle("Correlation with mu_1") +
  coord_flip() + theme_cowplot()
ggsave("figures/correlatedvariancemu1.png", height = 4.5, width = 8)


