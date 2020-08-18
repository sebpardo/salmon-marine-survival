# wrangling smolt data for Saint-Jean and Trinité rivers

library(tidyverse)
library(readxl)

source("analysis/99-hypergeometric-bootstrap-function.R")

# smolts data
# qcsmolts <- read_excel("raw-data/Smolt_St-Jean-Trinite.xlsx")
# names(qcsmolts) <- c("river", "year", "smolts")
# 
# table(qcsmolts$river)
# 
# ggplot(qcsmolts, aes(year, smolts, color = river)) + geom_point() + geom_line() + ylim(0,160000)
# 
# saveRDS(qcsmolts, file = "data/qc-smolts.rds")

# returns data
qcdata <- read_excel("raw-data/données_multifrayeur_Trinité et St-Jean_Quebec.xlsx")

tail(qcdata, n = 7)
qcdata <- qcdata[1:(nrow(qcdata) - 6), ] # removing last 6 rows
tail(qcdata, n = 7)

qcadults <- qcdata %>% select(Rivière, 
                     Année, 
                     `Montaison Mad`, 
                     `Montaison Réd`, 
                     `Montaison totale`,
                     `Age Nb échan. Diber`,
                     `Age Nb échan. Triber`,
                     `Age Nb échan. Fraie antérieure`,
                     `Age Nb en montaison Diber`, 
                     `Age Nb en montaison Triber`, 
                     `Age Nb en montaison Fraie antérieure`, 
                     `Reproducteur Mad Nb`, 
                     `Reproducteur Mad % Femelle`, 
                     `Reproducteur Mad Femelle Estimation`, 
                    `Reproducteur Réd Nb`, 
                    `Reproducteur Réd % Femelle`,	
                    `Reproducteur Réd Femelle Estimation`) %>%
          rename(river = Rivière, 
                 year = Année, 
                 returns_grilse = `Montaison Mad`, 
                 returns_msw = `Montaison Réd`, 
                 returns_total = `Montaison totale`,
                 scales_2sw = `Age Nb échan. Diber`,
                 scales_3sw = `Age Nb échan. Triber`,
                 scales_rs = `Age Nb échan. Fraie antérieure`,
                 returns_2sw_est = `Age Nb en montaison Diber`, 
                 returns_3sw_est = `Age Nb en montaison Triber`, 
                 returns_rs_est = `Age Nb en montaison Fraie antérieure`, 
                 spawners_grilse = `Reproducteur Mad Nb`, 
                 spawners_grilse_perc_female = `Reproducteur Mad % Femelle`, 
                 spawners_grilse_female_est = `Reproducteur Mad Femelle Estimation`,
                 spawners_msw = `Reproducteur Réd Nb`, 
                 spawners_msw_perc_female = `Reproducteur Réd % Femelle`,	
                 spawners_msw_female_est = `Reproducteur Réd Femelle Estimation`)

qcscales <- qcadults %>%
  filter(!(river == "Saint-Jean" & year == 1980)) %>% 
  select(river, year, returns_grilse, returns_2sw_est, returns_msw,
         scales_2sw, scales_3sw, scales_rs, returns_total) %>%
  mutate(n_small = returns_grilse, 
         n_large = returns_msw, 
         scales_small_1SW = returns_grilse, 
         scales_small_2SW = 0,
         scales_small_other = 0,
         scales_large_1SW = 0, 
         scales_large_2SW = scales_2sw, 
         scales_large_other = scales_3sw + scales_rs) %>%
  select(river,year,n_small:scales_large_other, returns_total) %>%
  mutate_all(replace_na, 0)

qcscalessds <-  pmap(qcscales, hyper_boot) %>% bind_rows %>%
  select(year, starts_with("sd")) %>%
  mutate(river_name = paste(qcscales$river, "River")) %>%
  select(river_name, everything())


# replacing NaNs (from years without samples) with mean value
qcscalessds[is.nan(qcscalessds$sdlog2SW),"sdlog2SW"] <- 0.999 #  mean(wabsds$sdlog2SW, na.rm = TRUE)

qcepsilon <- qcscales %>%
  mutate(river = paste(river, "River")) %>%
  select(river_name = river, year, n_small, n_large, 
                   totalreturns = returns_total, scales_small_1SW:scales_large_other) %>%
  left_join(qcscalessds, by = c("river_name","year"))

saveRDS(qcepsilon, file = "data/qcepsilon.rds")





saveRDS(filter(qcscalessds, river_name == "Trinité River"), "data/trinite-scales-sd.rds")
saveRDS(filter(qcscalessds, river_name == "Saint-Jean River"), "data/sj-scales-sd.rds")




qclong <- gather(qcadults, stage, numbers, -year, -river)


qclong %>%
  filter(grepl("^returns", stage)) %>%
  filter(stage %in% c("returns_grilse", "returns_2sw_est", "returns_3sw_est", "returns_rs_est", "returns_total")) %>%
ggplot( aes(year, numbers, color = stage)) + geom_point() + geom_line() + ylim(0, 3500) + facet_wrap(~river)

# Looking at estimated repeat spawners only
qclong %>%
  filter(stage == "returns_rs_est") %>%
  ggplot( aes(year, numbers, color = river)) + geom_point() + geom_line()

qcsmoltslong <- gather(qcsmolts, stage, numbers, -year, -river)

qcall <- bind_rows(qclong, qcsmoltslong)

# saving combined tidy dataset
saveRDS(qcall, file = "data/qc-stjean-trinite-all.rds")

