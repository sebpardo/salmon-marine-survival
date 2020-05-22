library(tidyverse)
library(rstan)
library(xtable)

fithierncp <- readRDS("data/stanfit_hier_ncp.rds") # hierarchical

alldatastan <- readRDS("data/alldatastan.rds")
river_names <- rep(alldatastan$river_name, times = alldatastan$nyears)

logsmoltsmat <- as.matrix(fithierncp, pars = "logsmolts_true")
logsmoltstib <- logsmoltsmat %>% t %>% as_tibble %>% mutate(river_name = river_names, year = alldatastan$years) %>%
  select(river_name, year, everything())

estsmolts <- pivot_longer(logsmoltstib, V1:V7500, names_to = "iteration") %>%
  mutate(year = year - 1) %>%
  group_by(river_name, year) %>%
  rename(`Population` = river_name, Year = year) %>%
  summarise(Median = exp(median(value)), 
            `2.5%` = exp(quantile(value, probs = 0.025)), 
            `25%` = exp(quantile(value, probs = 0.25)), 
            `75%` = exp(quantile(value, probs = 0.75)), 
            `97.5%` = exp(quantile(value, probs = 0.975))
           # n = n()
            )

allyears <- estsmolts %>%
  group_by(Population) %>%
  expand(Year = full_seq(Year, 1))

estsmolts_mis <- right_join(estsmolts, allyears, by = c("Population", "Year")) %>%
  mutate(Year = as.character(Year))

print.xtable(xtable(estsmolts_mis, digits = c(0,0,0,1,1,1,1,1),
                    display= c("s", "s", "s", rep("f", 5)),
                    caption = "Posterior annual smolt estimates for the seven populations examined. These
                    estimates differ from the smolt data used as they include both observation error and shrinkage from
                    the hierarchical model specification."),
            file = "ms/tableS1_smolts.tex",
            # only.contents = TRUE,
            NA.string = "-",
            tabular.environment = "longtable", floating = FALSE,
            size = "footnotesize",
            format.args = list(big.mark = ","),
            caption.placement = "top",
       include.rownames = FALSE)


logsmoltsmat <- as.matrix(fithierncp, pars = "logsmolts_mean")

mu_smolts <- as.matrix(fithierncp, pars = "logsmolts_mean") %>% 
  t %>% as_tibble %>% mutate(river_name = alldatastan$river_name) %>%
  pivot_longer(V1:V7500, names_to = "iteration") %>%
  select(river_name, everything()) %>%
  group_by(river_name) %>%
  rename(`Population` = river_name) %>%
  summarise(Median = exp(median(value)), 
            `2.5%` = exp(quantile(value, probs = 0.025)), 
            `25%` = exp(quantile(value, probs = 0.25)), 
            `75%` = exp(quantile(value, probs = 0.75)), 
            `97.5%` = exp(quantile(value, probs = 0.975))
            # n = n()
  )

mu_logsmolts <- as.matrix(fithierncp, pars = "logsmolts_mean") %>% 
  t %>% as_tibble %>% mutate(river_name = alldatastan$river_name) %>%
  pivot_longer(V1:V7500, names_to = "iteration") %>%
  select(river_name, everything()) %>%
  group_by(river_name) %>%
  rename(`Population` = river_name) %>%
  summarise(Median = (median(value)), 
            `2.5%` = (quantile(value, probs = 0.025)), 
            `25%` = (quantile(value, probs = 0.25)), 
            `75%` = (quantile(value, probs = 0.75)), 
            `97.5%` = (quantile(value, probs = 0.975))
            # n = n()
  )


sigma_logsmolts <- as.matrix(fithierncp, pars = "logsmolts_sigma") %>% 
  t %>% as_tibble %>% mutate(river_name = alldatastan$river_name) %>%
  pivot_longer(V1:V7500, names_to = "iteration") %>%
  select(river_name, everything()) %>%
  group_by(river_name) %>%
  rename(`Population` = river_name) %>%
  summarise(Median = median(value), 
            `2.5%` = quantile(value, probs = 0.025), 
            `25%` = quantile(value, probs = 0.25), 
            `75%` = quantile(value, probs = 0.75), 
            `97.5%` = quantile(value, probs = 0.975)
            # n = n()
  )

print.xtable(xtable(mu_smolts, digits = c(0,0,1,1,1,1,1),
                    caption = "Posterior estimates of the population-level 
                    mean smolt abundance for the seven populations examined."),
             file = "ms/tableS2_smu.tex",
             caption.placement = "top",
             format.args = list(big.mark = ","),
             include.rownames = FALSE)

print.xtable(xtable(mu_logsmolts, digits = c(0,0,3,3,3,3,3),
                    caption = "Posterior estimates of the population-level 
                    mean log-smolt abundance for the seven populations examined."),
             file = "ms/tableS3_lsmu.tex",
             caption.placement = "top",
             include.rownames = FALSE)

print.xtable(xtable(sigma_logsmolts, digits = c(0,0,3,3,3,3,3),
                    caption = "Posterior estimates of the population-level standard 
                    deviation of log-smolt abundance for the seven populations examined."),
             file = "ms/tableS4_lssig.tex",
             caption.placement = "top",
             include.rownames = FALSE)


# Tables of smolt estimates and cv values used


logsmolts <- tibble(river_name = rep(alldatastan$river_name, times = alldatastan$nyears),
                     year = alldatastan$years,
                     logsmolts = exp(alldatastan$logsmolts)
                   # logsmolts_cv = alldatastan$logsmolts_cv,
                   ) %>%
  arrange(river_name)

logsmolts_cv <- tibble(river_name = rep(alldatastan$river_name, times = alldatastan$nyears),
                    year = alldatastan$years,
                    #logsmolts = alldatastan$logsmolts,
                    logsmolts_cv = alldatastan$logsmolts_cv
) %>%
  arrange(river_name)

(smolts_tab <- pivot_wider(logsmolts, names_from = river_name, values_from = logsmolts) %>%
  arrange(year) %>%
  rename(`WAB` = "Western Arm Brook") %>%
  rename_all(. %>% gsub(" River", "", .)) %>%
  rename(Year = year) %>%
  mutate(Year = as.character(Year)))

print.xtable(xtable(smolts_tab, digits = c(0,0,rep(0,7)),
                    display= c("s", "s", rep("f",7)),
                    caption = "Empirical annual smolt estimates for the seven populations examined."),
             file = "ms/smolts_tab.tex",
             size = "footnotesize",
             format.args = list(big.mark = ","),
             tabular.environment = "longtable", floating = FALSE,
             NA.string = "-",
             caption.placement = "top",
             include.rownames = FALSE)



(smoltscv_tab <- pivot_wider(logsmolts_cv, names_from = river_name, values_from = logsmolts_cv) %>%
    arrange(year) %>%
    rename(`WAB` = "Western Arm Brook") %>%
    rename_all(. %>% gsub(" River", "", .)) %>%
    mutate_at(.vars = vars(-year), .funs = list(~ . * 100)) %>%
    rename(Year = year) %>%
    mutate(Year = as.character(Year)))

print.xtable(xtable(smoltscv_tab, digits = c(0,0,rep(1,7)),
                    caption = "Annual coefficient of variation (CV, \\%) estimates for the seven populations
                    examined. Where possible, we estimated CV directly
                    from uncertainty in smolts estimates."),
             file = "ms/smoltscv_tab.tex",
             NA.string = "-",
             size = "footnotesize",
             tabular.environment = "longtable", floating = FALSE,
             caption.placement = "top",
             include.rownames = FALSE)




