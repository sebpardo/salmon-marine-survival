### LaHave

library(tidyverse)
library(readxl)
library(lubridate)

# load biological characteristic data for returns (to assign 1SW-2SW)
lahave_bio_adult <- read_csv("raw-data/lahave/Lahave_TRadult.csv") %>%
  mutate(RCDATE = mdy(RCDATE))

# Reading additional 2018 data not needed because there aren't any smolt 
# abundance estimates for that year

lhbiochar <- lahave_bio_adult %>%
  select(-X1, -LOGID, -SEXID, -METHID, -CSITE, -GEARID, 
         -ORIGINID, -OTYPEID, -SCALE, -TAG, -EST, -BASIN) %>%
  mutate(RS = grepl("[Xx]", paste0(SP1, SP2, SP3, SP4, SP5))) %>%
  select(-SALMON) %>%
  filter(!is.na(POST)) %>% # removing those with no sea ages
  mutate(LS = ifelse(RS == TRUE, "RS", paste0(POST, "SW"))) 

# There are more large salmon scale samples than wild large retruns this year, so need to subset scales
lh1995 <- lahave_bio_adult %>%
  filter(RCYEAR == 1995)

# Test that only Xs are present
lhbiochar %>% 
  mutate(RS2 = paste0(SP1, SP2, SP3, SP4, SP5)) %>%
  pull(RS2) %>% 
  gsub("[NA]", "", .) %>% table

lhbiochar %>% filter(LS == "NASW")

# test that there are no false positives
lhbiochar %>% filter(RS == FALSE & !is.na(SP1))
lhbiochar %>% filter(RS == FALSE & !is.na(SP2))
lhbiochar %>% filter(RS == FALSE & !is.na(SP3))
lhbiochar %>% filter(RS == FALSE & !is.na(SP4))
lhbiochar %>% filter(RS == FALSE & !is.na(SP5))
lhbiochar %>% filter(RS == FALSE & !is.na(SP6))

lhbiowide <- lhbiochar %>% 
  select(-(SP1:SP6)) %>%
  rename(year = RCYEAR) %>%
  mutate(smla = ifelse(FORK >= 630, "large", "small")) %>%
  group_by(year, smla, LS) %>%
  summarise(number = n()) %>%
  spread(LS, number) %>%
  mutate(total = sum(`1SW`, `2SW`, `3SW`, `4SW`, `5SW`, RS, na.rm = TRUE)) %>%
  mutate(prop1SW = `1SW`/total, prop2SW = `2SW`/total, propRS = RS/total) %>%
  mutate(prop1SW = ifelse(is.na(prop1SW), 0, prop1SW)) %>%
  mutate(prop2SW = ifelse(is.na(prop2SW), 0, prop2SW)) %>%
  mutate(propRS = ifelse(is.na(propRS), 0, propRS)) %>%
  filter(year >= 1995)

# Reading return counts of WILD FISH ONLY!!!
lhreturns <- read_excel("raw-data/lahave/Morgan Falls adult returns.xlsx", 
                        skip = 2, na = c("", "--"), n_max = 50) %>%
  select(year = `...1`, small = `Wild Small`, large = `Wild Large`)

lhreturnslong <- lhreturns %>%
  gather(key = smla, value = returns, -year) %>%
  arrange(year, desc(smla)) %>%
  filter(year >= 1995)

lhreturns_corrected <- left_join(lhreturnslong, lhbiowide, by = c("year", "smla")) %>%
  mutate(est1SW =  returns * prop1SW) %>%
  mutate(est2SW = returns * prop2SW) %>%
  mutate(estRS = returns * propRS) %>%
  select(year, smla, returns, est1SW, est2SW, estRS) %>%
  ungroup() %>%
  group_by(year) %>%
  arrange(year, desc(smla)) %>%
  summarise(est1SW = sum(est1SW), 
            est2SW = sum(est2SW), 
            estRS = sum(estRS),
            totalsmall = first(returns),
            totallarge = last(returns))



# Group numbers of scales for CV calculation
lhwide <- lhbiowide %>%
  select(year, smla, `1SW`:RS) %>% 
  pivot_wider(names_from = smla, values_from = c(`1SW`:RS))

lhwide[is.na(lhwide)] <- 0 # replacing NAs with zeros

lhjoin <- left_join(lhreturns_corrected, lhwide, by = "year")

lhscales <- lhjoin %>%
  mutate_all(replace_na, 0) %>%
  mutate(n_small = totalsmall, 
         n_large = totallarge,
         scales_small_1SW = `1SW_small`,
         scales_small_2SW = `2SW_small`,
         scales_small_other = `3SW_small` + `4SW_small` + `5SW_small` + RS_small,
         scales_large_1SW = `1SW_large`, 
         scales_large_2SW = `2SW_large`, 
         scales_large_other = `3SW_large` + `4SW_large` + `5SW_large` + RS_large) %>%
  select(year, n_small:scales_large_other)

# NOT WORKING, more scale samples that wild returns in a few years (starting 1995) suggesting
# that scales come from both wild and hatchery returns
# lhsds <-  pmap(lhscales, hyper_boot) %>% bind_rows %>%
#   select(year, starts_with("sd"))


# Reading smolt abundance data
lhsmolts <- read_excel("raw-data/lahave/Morgan Falls Smolt estimates.xlsx") %>%
  rename(year = Year, smolts_med = Estimate, smolts_lower = Lower, smolts_upper = Upper)

# Reading smolt biological characteristics data
lhsmolts_biochar <- read_excel("raw-data/lahave/Morgan Falls Smolt Biol Characteristics.xlsx",
                               na = c("", "NA"),
                               col_types = c(rep("guess", 7), "text", "guess", "guess", rep("text", 3))) %>%
  rename(ra = `Fw.Age`) %>%
  select(Year:ra) %>%
  mutate(date = as_date(paste(Year, Month, Day, sep = "-")))


# Merging smolt and return estimates
lahavedat <- left_join(lhreturns_corrected, lhsmolts, by = "year")


lahavereturns <- lahavedat %>%
  mutate(smoltdiff =  ((smolts_med-smolts_lower)+(smolts_upper-smolts_med))/2,
         smoltdiff_upper =  smolts_upper-smolts_med,
         smoltdiff_lower =  smolts_med-smolts_lower,
         smoltsd = smoltdiff/1.96,
         smoltscv = smoltsd/smolts_med,
         smoltlag = lag(smolts_med),
         smolt025lag = lag(smolts_lower),
         smolt975lag = lag(smolts_upper),
         est_1SW = est1SW,
         est_2SWlead = lead(est2SW),
         est_RSlead = lead(estRS),
         smoltscvlag = lag(smoltscv)) 

# which years have missing smolt data?

# lahavereturns[!is.na(lahavereturns$smoltlag),]

# removing years missing smolt data (lagged)
lahavereturns <- lahavereturns %>%
  filter(!is.na(smoltlag))


# Creating list for use in TMB
lahavedata <- list(years = lahavereturns$year,
                     logsmolts = log(lahavereturns$smoltlag),
                     logsmolts_cv = lahavereturns$smoltscvlag,
                     loggrilse = log(lahavereturns$est_1SW),
                     logSW2 = log(lahavereturns$est_2SWlead),
                     logRS = log(lahavereturns$est_RSlead),
                     #N =  length(connereturns$small[-1]),
                     #Pr = rep(0.8, nrow(connereturns)))
                     returns_cv = 0.01, # ASSUMPTION, CV of return abundance counts
                   logSW1_cv = rep(0.01, length(lahavereturns$year)), # NO DATA SO USING VALUES BASED ON MEAN TRINITE ESTIMATES
                   logSW2_cv =  rep(0.06, length(lahavereturns$year)),
                     nyears = length(lahavereturns$year),
                     river_name = "LaHave River",
                     allreturns = lahavedat %>% filter(year >= 1995) %>% mutate(river_name = "LaHave River")
) 

saveRDS(lahavedata, file = "data/lahavedata.rds")


