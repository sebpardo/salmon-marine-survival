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

# Writing csv of (unlogged) smolts_true for use as input in CPM
write_csv(estsmolts, file = "data/smolts-true-posteriors.csv")

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
    mutate(year = year - 1) %>%
    rename(`WAB` = "Western Arm Brook") %>%
  rename_all(. %>% gsub(" River", "", .)) %>%
  rename(Year = year) %>%
  mutate(Year = as.character(Year)))

print.xtable(xtable(smolts_tab, digits = c(0,0,rep(0,7)),
                    display= c("s", "s", rep("f",7)),
                    caption = "Empirically estimated annual smolt abundances for the seven populations examined. Smolt abundance estimates from the Trinité, Saint-Jean, LaHave, Nashwaak, and Conne populations were obtained using a mark-recapture approach, while estimates from the WAB and Campbellton populations were obtained by direct counts using fish counting fences."),
             file = "ms/smolts_tab.tex",
             size = "footnotesize",
             format.args = list(big.mark = ","),
             tabular.environment = "longtable", floating = FALSE,
             NA.string = "-",
             caption.placement = "top",
             include.rownames = FALSE)



(smoltscv_tab <- pivot_wider(logsmolts_cv, names_from = river_name, values_from = logsmolts_cv) %>%
    arrange(year) %>%
    mutate(year = year - 1) %>%
    rename(`WAB` = "Western Arm Brook") %>%
    rename_all(. %>% gsub(" River", "", .)) %>%
    mutate_at(.vars = vars(-year), .funs = list(~ . * 100)) %>%
    rename(Year = year) %>%
    mutate(Year = as.character(Year)))

print.xtable(xtable(smoltscv_tab, digits = c(0,0,rep(1,7)),
                    caption = "Annual coefficient of variation (CV, \\%) in smolt estimates for the seven populations
                    examined. Where possible, we estimated CV directly
                    from uncertainty in smolts estimates."),
             file = "ms/smoltscv_tab.tex",
             NA.string = "-",
             size = "footnotesize",
             tabular.environment = "longtable", floating = FALSE,
             caption.placement = "top",
             include.rownames = FALSE)


##### EPSILON ESTIMATES TABLE

conneepsilon <- readRDS("data/conneepsilon.rds") %>% filter(year %in% na.omit(alldatastan$yearmat[3,]))
wabepsilon <- readRDS("data/wabepsilon.rds")  %>% filter(year %in% na.omit(alldatastan$yearmat[5,]))
sjepsilon <- readRDS("data/qcepsilon.rds") %>% 
  filter(river_name == "Saint-Jean River" & year %in% na.omit(alldatastan$yearmat[1,]))
triepsilon <- readRDS("data/qcepsilon.rds") %>% 
  filter(river_name == "Trinité River" & year %in% na.omit(alldatastan$yearmat[2,]))
campepsilon <- readRDS("data/campepsilon.rds") %>% filter(year %in% na.omit(alldatastan$yearmat[6,]))

allscales <- bind_rows(campepsilon, wabepsilon, conneepsilon, sjepsilon, triepsilon)

epsilonslh <- tibble(river_name = "LaHave River",  
                     year = na.omit(alldatastan$yearmat[7,]),
                     sdlog1SW = 0.01, sdlog2SW = 0.06)

epsilonsnash <- tibble(river_name = "Nashwaak River",  
                     year = na.omit(alldatastan$yearmat[4,]),
                     sdlog1SW = 0.01, sdlog2SW = 0.06)


epsilons <- bind_rows(allscales, epsilonslh, epsilonsnash) %>%
  select(river_name, year, sdlog1SW, sdlog2SW)

epsilons1 <- select(epsilons, river_name, year, sdlog1SW)
epsilons2 <- select(epsilons, river_name, year, sdlog2SW)



(epsilon1_tab <- epsilons1 %>%
    arrange(river_name) %>% 
    pivot_wider(names_from = river_name, values_from = sdlog1SW) %>%
    arrange(year) %>%
    rename(`WAB` = "Western Arm Brook") %>%
    rename_all(. %>% gsub(" River", "", .)) %>%
    # mutate_at(.vars = vars(-year), .funs = list(~ . * 100)) %>%
    rename(Year = year) %>%
    mutate(Year = as.character(Year)))

print.xtable(xtable(epsilon1_tab, digits = c(0,0,rep(3,7)),
                    caption = "Annual error estimates derived from the standard deviation of log-transformed abundance 
                    in 1SW returns ($\\epsilon_{1}$), derived from size group to age conversion using aged scale data 
                    (except for Nashwaak and LaHave populations).",
                    label = "tab:epsilon1"),
             file = "ms/epsilon1_tab.tex",
             NA.string = "-",
             size = "footnotesize",
             tabular.environment = "longtable", floating = FALSE,
             caption.placement = "top",
             include.rownames = FALSE)


(epsilon2_tab <- epsilons2 %>%
    mutate(sdlog2SW = if_else(sdlog2SW == 0.99, 0.999, sdlog2SW)) %>%
    arrange(river_name) %>% 
    pivot_wider(names_from = river_name, values_from = sdlog2SW) %>%
    arrange(year) %>%
    mutate(year = year + 1) %>%
    rename(`WAB` = "Western Arm Brook") %>%
    rename_all(. %>% gsub(" River", "", .)) %>%
    #mutate_at(.vars = vars(-year), .funs = list(~ . * 100)) %>%
    rename(Year = year) %>%
    mutate(Year = as.character(Year)))





print.xtable(xtable(epsilon2_tab, digits = c(0,0,rep(3,7)),
                    caption = "Annual error estimates derived from the standard deviation of log-transformed abundance 
                    in 2SW returns ($\\epsilon_{2}$), derived from size group to age conversion using aged scale data 
                    (except for Nashwaak and LaHave populations).",
                    label = "tab:epsilon2"),
             file = "ms/epsilon2_tab.tex",
             NA.string = "-",
             size = "footnotesize",
             tabular.environment = "longtable", floating = FALSE,
             caption.placement = "top",
             include.rownames = FALSE)


allscales_tab <- allscales %>%
  select(Population = river_name, 
         Year = year, 
         #'Small returns' = n_small, 
         #'Large returns' = n_large, 
         #'Total returns' = totalreturns,
         'Scales small 1SW' = scales_small_1SW,
         'Scales small 2SW' = scales_small_2SW,
         'Scales small other' = scales_small_other,
         'Scales large 1SW' = scales_large_1SW,
         'Scales large 2SW' = scales_large_2SW,
         'Scales large other' = scales_large_other)
         
print.xtable(xtable(allscales_tab, digits = rep(0,9),
                    align = c("p{0cm}","p{3cm}","p{1cm}",rep("p{1.3cm}",6)),
                    caption = "Scale data used to estimate $\\epsilon_{1}$ 
                    and $\\epsilon_{2}$ for five of the seven populations examined (we could not obtain
                    scale data for Nashwaak and LaHave populations). The ``other'' category includes
                    maiden spawners of sea age 3 or more and repeat spawners.",
                    label = "tab:scales"),
             file = "ms/scales_tab.tex",
             NA.string = "-",
             size = "footnotesize",
             tabular.environment = "longtable", floating = FALSE,
             caption.placement = "top",
             include.rownames = FALSE)
