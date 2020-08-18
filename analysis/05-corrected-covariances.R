# Corrected covariances between mu1 and estimated parameters 

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


# mu1var %>%
#   mutate(river_name = fct_relevel(river_name, propgrilse$river_name),
#          correlation = fct_relevel(correlation, "mu sum", after = Inf)) %>%
#   ggplot(aes(correlation, median, color = river_name)) + 
#   geom_hline(yintercept = 0, col = "gray80") +
#   geom_point(position = position_dodge(-0.5)) + 
#   geom_errorbar(aes(ymax = upperCI, ymin = lowerCI), width = 0.0, position = position_dodge(-0.5)) +
#   ylab("Contribution to Variance") + xlab("Parameter") + ggtitle("Correlation with mu_1") +
#   coord_flip() + theme_cowplot()
# ggsave("figures/correlatedvariancemu1.png", height = 4.5, width = 8)



mu2var <- tibble(river_name = rep(unique(river_names), each = 11),
                 correlation = rep(c("Var(logsmolts)/Var(mu)", 
                                     "Var(Z1)/Var(mu)", 
                                     "Var(Z2)/Var(mu)", 
                                     "Var(logPr)/Var(mu)",
                                     "2*Cov(logsmolts,-Z1)/Var(mu)",
                                     "2*Cov(logsmolts,-Z2)/Var(mu)",
                                     "2*Cov(logsmolts,logPr)/Var(mu)",
                                     "2*Cov(-Z1,logPr)/Var(mu)",
                                     "2*Cov(-Z2,logPr)/Var(mu)",
                                     "2*Cov(-Z1,-Z2)/Var(mu)",
                                     "mu sum"), 7),
                 median = NaN,
                 lowerCI = NaN,
                 upperCI = NaN)

mu2errorvar <- tibble(river_name = rep(unique(river_names), each = 11),
                      correlation = rep(c("Var(logsmolts)", 
                                          "Var(Z1)", 
                                          "Var(Z2)", 
                                          "Var(logPr)",
                                          "Var(mu)",
                                          "Cov(logsmolts,-Z1)",
                                          "Cov(logsmolts,-Z2)",
                                          "Cov(logsmolts,logPr)",
                                          "Cov(-Z1,logPr)",
                                          "Cov(-Z2,logPr)",
                                          "Cov(-Z1,-Z2)"), 7),
                      mean = NaN)

mu2correctedcov <- tibble(river_name = rep(unique(river_names), each = 6),
                          correlation = rep(c("Corrected cov(logsmolts,-Z1)",
                                              "Corrected cov(logsmolts,-Z2)",
                                              "Corrected cov(logsmolts,logPr)",
                                              "Corrected cov(-Z1,logPr)",
                                              "Corrected cov(-Z2,logPr)",
                                              "Corrected cov(-Z1,-Z2)"), 7),
                          median = NaN,
                          lowerCI = NaN,
                          upperCI = NaN)

mu2correctedrsq <- tibble(river_name = rep(unique(river_names), each = 5),
                          correlation = rep(c("R2(mu2,-Z1)",
                                              "R2(mu2,-Z2)",
                                              "R2(mu2,logsmolts)",
                                              "R2(mu2,logPr)",
                                              "R2 sum"), 7),
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
  varlogpr <- mapply(var, log(prsub))
  varz1 <- mapply(var, z1sub)
  varz2 <- mapply(var, z2sub)
  varlogsmolts <- mapply(var, logsmoltssub)
  #browser()
  
  covlogsmoltsz1 <- mapply(cov, logsmoltssub, z1sub*(-1)) 
  covlogsmoltsz2 <- mapply(cov, logsmoltssub, z2sub*(-1)) 
  covlogsmoltslogpr <- mapply(cov, logsmoltssub, log(prsub)) 
  covz1logpr <- mapply(cov, z1sub*(-1),  log(prsub)) 
  covz2logpr <- mapply(cov, z2sub*(-1),  log(prsub)) 
  covz1z2 <- mapply(cov, z1sub*(-1),  z2sub*(-1)) 
  
  covmu2logsmolts <- mapply(cov, mu2sub, logsmoltssub) 
  covmu2z1 <- mapply(cov, mu2sub, z1sub*(-1)) 
  covmu2z2 <- mapply(cov, mu2sub, z2sub*(-1)) 
  covmu2logpr <- mapply(cov, mu2sub,  log(prsub)) 
  
    # estimating error variances and covariances
  # this is done by estimating variance across iterations for each year and then averaging across years,
  # hence the transpose
  errorvarmu2 <- mapply(var, mu2sub %>% t %>% as_tibble) %>% mean
  errorvarlogpr <- mapply(var, log(prsub) %>% t %>% as_tibble) %>% mean
  errorvarz1 <- mapply(var, z1sub %>% t %>% as_tibble) %>% mean
  errorvarz2 <- mapply(var, z2sub %>% t %>% as_tibble) %>% mean
  errorvarlogsmolts <- mapply(var, logsmoltssub %>% t %>% as_tibble) %>% mean
  
  errorcovlogsmoltsz1 <- mapply(cov, logsmoltssub %>% t %>% as_tibble, t(z1sub*(-1)) %>% as_tibble)  %>% mean
  errorcovlogsmoltsz2 <- mapply(cov, logsmoltssub %>% t %>% as_tibble, t(z2sub*(-1)) %>% as_tibble)  %>% mean
  errorcovlogsmoltslogpr <- mapply(cov, logsmoltssub %>% t %>% as_tibble, log(prsub) %>% t %>% as_tibble)  %>% mean
  errorcovz1logpr <- mapply(cov, t(z1sub*(-1)) %>% as_tibble,  log(prsub) %>% t %>% as_tibble)  %>% mean
  errorcovz2logpr <- mapply(cov, t(z2sub*(-1)) %>% as_tibble,  log(prsub) %>% t %>% as_tibble)  %>% mean
  errorcovz1z2 <- mapply(cov, t(z1sub*(-1)) %>% as_tibble,  t(z2sub*(-1)) %>% as_tibble)  %>% mean
  
  errorcovmu2logsmolts <- mapply(cov, mu2sub %>% t %>% as_tibble,logsmoltssub %>% t %>% as_tibble)  %>% mean
  errorcovmu2z1 <- mapply(cov, mu2sub %>% t %>% as_tibble, t(z1sub*(-1)) %>% as_tibble)  %>% mean
  errorcovmu2z1 <- mapply(cov, mu2sub %>% t %>% as_tibble, t(z2sub*(-1)) %>% as_tibble)  %>% mean
  errorcovmu2logpr <- mapply(cov, mu2sub %>% t %>% as_tibble, log(prsub) %>% t %>% as_tibble)  %>% mean
  
  # saving to out
  mu2errorvar[mu2errorvar$river_name == i & mu2errorvar$correlation == "Var(mu)", 3] <- errorvarmu2
  mu2errorvar[mu2errorvar$river_name == i & mu2errorvar$correlation == "Var(logPr)", 3] <- errorvarlogpr
  mu2errorvar[mu2errorvar$river_name == i & mu2errorvar$correlation == "Var(Z1)", 3] <- errorvarz1
  mu2errorvar[mu2errorvar$river_name == i & mu2errorvar$correlation == "Var(Z2)", 3] <- errorvarz2
  mu2errorvar[mu2errorvar$river_name == i & mu2errorvar$correlation == "Var(logsmolts)", 3] <- errorvarlogsmolts
  
  mu2errorvar[mu2errorvar$river_name == i & mu2errorvar$correlation == "Cov(logsmolts,-Z1)", 3] <- errorcovlogsmoltsz1
  mu2errorvar[mu2errorvar$river_name == i & mu2errorvar$correlation == "Cov(logsmolts,-Z2)", 3] <- errorcovlogsmoltsz2
  mu2errorvar[mu2errorvar$river_name == i & mu2errorvar$correlation == "Cov(logsmolts,logPr)", 3] <- errorcovlogsmoltslogpr
  mu2errorvar[mu2errorvar$river_name == i & mu2errorvar$correlation == "Cov(-Z1,logPr)", 3] <- errorcovz1logpr
  mu2errorvar[mu2errorvar$river_name == i & mu2errorvar$correlation == "Cov(-Z2,logPr)", 3] <- errorcovz2logpr
  mu2errorvar[mu2errorvar$river_name == i & mu2errorvar$correlation == "Cov(-Z1,-Z2)", 3] <- errorcovz1z2
  
  # mu2errorvar[mu2errorvar$river_name == i & mu2errorvar$correlation == "Cov(mu1,logsmolts)", 3] <- errorcovmu1logsmolts
  # mu2errorvar[mu2errorvar$river_name == i & mu2errorvar$correlation == "Cov(mu1,-Z1)", 3] <- errorcovmu1z1
  # mu2errorvar[mu2errorvar$river_name == i & mu2errorvar$correlation == "Cov(mu1,-Z2)", 3] <- errorcovmu1z2
  # mu2errorvar[mu2errorvar$river_name == i & mu2errorvar$correlation == "Cov(mu1,logPr)", 3] <- errorcovmu1logpr
  
  # estimating corrected correlation values
  corrected_corlogsmoltsz1 <- (covlogsmoltsz1 - errorcovlogsmoltsz1)/(sqrt(varlogsmolts - errorvarlogsmolts)*sqrt(varz1 - errorvarz1))
  corrected_corlogsmoltsz2 <- (covlogsmoltsz2 - errorcovlogsmoltsz1)/(sqrt(varlogsmolts - errorvarlogsmolts)*sqrt(varz2 - errorvarz2))
  corrected_corlogsmoltslogpr <- (covlogsmoltslogpr - errorcovlogsmoltslogpr)/(sqrt(varlogsmolts - errorvarlogsmolts)*sqrt(varlogpr - errorvarlogpr))
  corrected_corz1logpr <- (covz1logpr - errorcovz1logpr)/(sqrt(varlogpr - errorvarlogpr)*sqrt(varz1 - errorvarz1))
  corrected_corz2logpr <- (covz2logpr - errorcovz2logpr)/(sqrt(varlogpr - errorvarlogpr)*sqrt(varz2 - errorvarz2))
  corrected_corz1z2 <- (covz1z2 - errorcovz1z2)/(sqrt(varz1 - errorvarz1)*sqrt(varz2 - errorvarz2))
  
  # estimating R^2 with mu1
  # corrected_r2_mu1logsmolts <- ((covmu1logsmolts - errorcovmu1logsmolts)/(sqrt(varlogsmolts - errorvarlogsmolts)*sqrt(varmu1 - errorvarmu1)))^2
  # corrected_r2_mu1z1 <- ((covmu1z1 - errorcovmu1z1)/(sqrt(varz1 - errorvarz1)*sqrt(varmu1 - errorvarmu1)))^2
  # corrected_r2_mu1logpr <- ((covmu1logpr - errorcovmu1logpr)/(sqrt(varlogpr - errorvarlogpr)*sqrt(varmu1 - errorvarmu1)))^2
  # 
  # rsqsum <- colSums(rbind(corrected_r2_mu1logsmolts, corrected_r2_mu1z1, corrected_r2_mu1logpr))
  
  
  # storing values, na.rm = TRUE to ignore NaNs created when square-rooting a few negative values
  mu2correctedcov[mu2correctedcov$river_name == i & 
                    mu2correctedcov$correlation == "Corrected cov(logsmolts,-Z1)", 3:5] <- quantile(corrected_corlogsmoltsz1, 
                                                                                                    c(0.5, 0.05, 0.95), na.rm = TRUE)
  mu2correctedcov[mu2correctedcov$river_name == i & 
                    mu2correctedcov$correlation == "Corrected cov(logsmolts,-Z2)", 3:5] <- quantile(corrected_corlogsmoltsz2, 
                                                                                                    c(0.5, 0.05, 0.95), na.rm = TRUE)
  mu2correctedcov[mu2correctedcov$river_name == i & 
                    mu2correctedcov$correlation == "Corrected cov(logsmolts,logPr)", 3:5] <- quantile(corrected_corlogsmoltslogpr, 
                                                                                                      c(0.5, 0.05, 0.95), na.rm = TRUE)
  mu2correctedcov[mu2correctedcov$river_name == i & 
                    mu2correctedcov$correlation == "Corrected cov(-Z1,logPr)", 3:5] <- quantile(corrected_corz1logpr, 
                                                                                                c(0.5, 0.05, 0.95), na.rm = TRUE)
  mu2correctedcov[mu2correctedcov$river_name == i & 
                    mu2correctedcov$correlation == "Corrected cov(-Z2,logPr)", 3:5] <- quantile(corrected_corz2logpr, 
                                                                                                c(0.5, 0.05, 0.95), na.rm = TRUE)
  mu2correctedcov[mu2correctedcov$river_name == i & 
                    mu2correctedcov$correlation == "Corrected cov(-Z1,-Z2)", 3:5] <- quantile(corrected_corz1z2, 
                                                                                                c(0.5, 0.05, 0.95), na.rm = TRUE)
  
  # storing values, na.rm = TRUE to ignore NaNs created when square-rooting a few negative values
  # mu1correctedrsq[mu1correctedrsq$river_name == i & 
  #                   mu1correctedrsq$correlation == "R2(mu1,logsmolts)", 3:5] <- quantile(corrected_r2_mu1logsmolts,
  #                                                                                        c(0.5, 0.05, 0.95), na.rm = TRUE)
  # mu1correctedrsq[mu1correctedrsq$river_name == i & 
  #                   mu1correctedrsq$correlation == "R2(mu1,-Z1)", 3:5] <- quantile(corrected_r2_mu1z1,
  #                                                                                  c(0.5, 0.05, 0.95), na.rm = TRUE)
  # mu1correctedrsq[mu1correctedrsq$river_name == i & 
  #                   mu1correctedrsq$correlation == "R2(mu1,logPr)", 3:5] <- quantile(corrected_r2_mu1logpr,
  #                                                                                    c(0.5, 0.05, 0.95), na.rm = TRUE)
  # mu1correctedrsq[mu1correctedrsq$river_name == i & 
  #                   mu1correctedrsq$correlation == "R2 sum", 3:5] <- quantile(rsqsum, c(0.5, 0.05, 0.95), na.rm = TRUE)
  
  # varsmoltsmu <- varlogsmolts/varmu1
  # varz1mu <- varz1/varmu1
  # varlogprmu <- varlogpr/varmu1
  # cov2logsmoltsz1mu <- 2*covlogsmoltsz1/varmu1
  # cov2logsmoltslogprmu <- 2*covlogsmoltslogpr/varmu1
  # cov2z1logprmu <- 2*covz1logpr/varmu1
  # 
  # rsqsum <- colSums(rbind(varsmoltsmu, varz1mu, varlogprmu, cov2logsmoltsz1mu, cov2logsmoltslogprmu, cov2z1logprmu))
  
  print(i)
  # print(paste("R squared mu with Z1 =", median(mu1z1)))
  # print(paste("R squared mu with log(Pr) =", median(mu1pr)))
  # print(paste("R squared mu with logsmolts_true =", median(mu1logsmolts)))
  # print(paste("median R squared sum =", median(rsqsum)))
  
  # mu1var[mu1var$river_name == i & mu1var$correlation == "Var(logsmolts)/Var(mu)", 3] <- median(varsmoltsmu)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "Var(logsmolts)/Var(mu)", 4] <- quantile(varsmoltsmu, 0.05)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "Var(logsmolts)/Var(mu)", 5] <- quantile(varsmoltsmu, 0.95)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "Var(Z1)/Var(mu)", 3] <- median(varz1mu)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "Var(Z1)/Var(mu)", 4] <- quantile(varz1mu, 0.05)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "Var(Z1)/Var(mu)", 5] <- quantile(varz1mu, 0.95)  
  # mu1var[mu1var$river_name == i & mu1var$correlation == "Var(logPr)/Var(mu)", 3] <- median(varlogprmu)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "Var(logPr)/Var(mu)", 4] <- quantile(varlogprmu, 0.05)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "Var(logPr)/Var(mu)", 5] <- quantile(varlogprmu, 0.95)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(logsmolts,-Z1)/Var(mu)", 3] <- median(cov2logsmoltsz1mu)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(logsmolts,-Z1)/Var(mu)", 4] <- quantile(cov2logsmoltsz1mu, 0.05)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(logsmolts,-Z1)/Var(mu)", 5] <- quantile(cov2logsmoltsz1mu, 0.95)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(logsmolts,logPr)/Var(mu)", 3] <- median(cov2logsmoltslogprmu)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(logsmolts,logPr)/Var(mu)", 4] <- quantile(cov2logsmoltslogprmu, 0.05)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(logsmolts,logPr)/Var(mu)", 5] <- quantile(cov2logsmoltslogprmu, 0.95)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(-Z1,logPr)/Var(mu)", 3] <- median(cov2z1logprmu)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(-Z1,logPr)/Var(mu)", 4] <- quantile(cov2z1logprmu, 0.05)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "2*Cov(-Z1,logPr)/Var(mu)", 5] <- quantile(cov2z1logprmu, 0.95)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "mu sum", 3] <- median(rsqsum)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "mu sum", 4] <- quantile(rsqsum, 0.05)
  # mu1var[mu1var$river_name == i & mu1var$correlation == "mu sum", 5] <- quantile(rsqsum, 0.95)  
}

#saveRDS(mu2correctedcov, file = "data/corrected-correlations.rds")
#mu2correctedcov <- readRDS("data/corrected-correlations.rds")
library(scales)

mu2correctedcov

pivot_wider(mu2errorvar, names_from = correlation, values_from = mean) %>%
  select(river_name, starts_with("Cov"))

mu2correctedcov %>% pull(correlation) %>% unique

mu2correctedcov %>%
  select(-lowerCI, -upperCI) %>%
  mutate(correlation = str_replace(correlation, "cov", "cor")) %>%
  pivot_wider(names_from = correlation, values_from = median)

# mu2correctedrsq %>%
#   #mutate(correlation = str_replace(correlation, "cov", "cor")) %>%
#   mutate(river_name = fct_relevel(river_name, propgrilse$river_name)) %>%
#   ggplot(aes(correlation, median, color = river_name)) + 
#   geom_hline(yintercept = 1, col = "gray80") +
#   geom_point(position = position_dodge(-0.5)) + 
#   geom_errorbar(aes(ymax = upperCI, ymin = lowerCI), width = 0.0, position = position_dodge(-0.5)) +
#   ylab("R squared value") + xlab("Parameters") + ggtitle("Corrected R squares") +
#   labs(color = "River name") +
#   coord_flip() + theme_cowplot()
# ggsave("figures/corrected_rsq_mu1.png", height = 4.5, width = 8)


pcorrectedcor <- mu2correctedcov %>%
  mutate(lowerCI = squish(lowerCI, c(-1,1)),
         upperCI = squish(upperCI, c(-1,1))) %>%
  mutate(correlation = str_replace(correlation, "cov", "cor")) %>%
  mutate(correlation2 = fct_recode(correlation,
    "italic(cor(log(smolts),-Z[1]))" = "Corrected cor(logsmolts,-Z1)" ,
    "italic(cor(log(smolts),-Z[2]))" = "Corrected cor(logsmolts,-Z2)",
    "italic(cor(log(smolts),log(P[g])))" = "Corrected cor(logsmolts,logPr)",
    "italic(cor(-Z[1],log(P[g])))" = "Corrected cor(-Z1,logPr)",
    "italic(cor(-Z[2],log(P[g])))" = "Corrected cor(-Z2,logPr)",
    "italic(cor(-Z[1],-Z[2]))" = "Corrected cor(-Z1,-Z2)" 
  )) %>%
  mutate(river_name = fct_relevel(river_name, propgrilse$river_name)) %>%
  ggplot(aes(correlation2, median, color = river_name)) + 
  geom_hline(yintercept = 0, col = "gray80") +
  geom_hline(yintercept = -1, col = "gray80") +
  geom_hline(yintercept = 1, col = "gray80") +
  geom_point(position = position_dodge(-0.5)) + 
  #coord_cartesian(clip = 'off') +
  geom_errorbar(aes(ymax = upperCI, ymin = lowerCI), width = 0.0, position = position_dodge(-0.5)) +
  ylab("Correlation value") + xlab("Parameters") + 
  ggtitle("Corrected correlations between estimated parameters") +
  labs(color = "Population") +
  scale_x_discrete(labels = label_parse()) +
  scale_y_continuous(limits = c(-1,1), oob = function(x, ...) x) +
  coord_cartesian(clip = "off") +
  #ylim(-1, 1) +
  coord_flip() + theme_cowplot() 
pcorrectedcor
ggsave("figures/corrected-cor-ALL.png", height = 4.5, width = 8.5)

pcorrectedcor + 
  scale_y_continuous(breaks = seq(-5, 2, by = 1), oob = function(x, ...) x)
ggsave("figures/corrected-cor-ALL-fullerrorbars.png", height = 4.5, width = 8.5)
  

