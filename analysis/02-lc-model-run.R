library(tidyverse)
library(TMB)

# Compiling TMB model

#compile("lc_model.cpp", "-O0 -g") # verbose output for gdbsource
compile("lc_model.cpp")

# dynamic link to object
dyn.load(dynlib("lc_model"))

connereturnsdat <- readRDS(file = "~/Dropbox/salmon-bayesian-survival/connereturns.rds")

# Saving filtered data frame in new object (connereturns)
connereturns <- connereturnsdat %>%
  mutate(smoltlag = lag(smolt),
         largelead = lead(large)) %>% 
  filter(year %in% 1990:2015) 

# Creating list for use in TMB
connedata <- list(#years = connereturns$year
  logsmolts = log(connereturns$smoltlag),
  logsmolts_SE =  rep(0.05, length(connereturns$smoltlag)),
  loggrilse = log(connereturns$small),
  logSW2 = log(connereturns$largelead),
  #N =  length(connereturns$small[-1]),
  #Pr = rep(0.8, nrow(connereturns)))
  cv = 0.01 # ASSUMPTION
  ) 


# salmondata <- connedata
salmondata <- nashwaakdata

pars <- list(Z1 = rep(-log(0.08), length(salmondata$logsmolts)), 
             Z2 = rep(-log(0.1), length(salmondata$logsmolts)),
             logitPr =  rep(0, length(salmondata$logsmolts)),
             logsmolts_true =  rep(mean(salmondata$logsmolts), length(salmondata$logsmolts)),
             logsmolts_mean = mean(salmondata$logsmolts),
             logsmolts_logsd = sd(salmondata$logsmolts)/1.3,
             Z1_mean = 2.5,
             Z1_logsd = 1,
             Z2_mean = 2.5,
             Z2_logsd = 1,
             logitPr_mean = 0,
             logitPr_logsd = 1
             )
             #mu = rep(mean(salmondata$loggrilse), length(salmondata$smolts)),
             #mu2 =  rep(mean(salmondata$logSW2), length(salmondata$smolts)))
               
# checking the lengths of parameters, they should match expectations from TMB model
lapply(pars, length)
# gdbsource("lc_model_run.R")

# running model
obj <- MakeADFun(salmondata, pars, DLL = "lc_model", random = c("logsmolts_true", "Z1", "Z2", "logitPr"))

str(obj)
# optimizing model objective
res <- nlminb(obj$par, obj$fn, obj$gr)
res$objective
res$par

# Calculating SDs using Laplace approximation
results <- sdreport(obj)

summary(results, select = "fixed")
estsd <- summary(results, select = "random")
estsd <- tibble::rownames_to_column(as.data.frame(estsd))
estsd$year <- rep(salmondata$year, nrow(estsd)/length(salmondata$year))

estsd <- estsd %>% as_tibble() %>%
  rename(parameter = rowname,
         value = Estimate,
         stderr = `Std. Error`) %>%
  select(parameter, year, everything())


estsd %>%
  #filter(grepl("^Z", parameter)) %>%
  mutate(lower = value - stderr,
         upper = value + stderr,
         value = ifelse(grepl("^Z", parameter), exp(-value), value),
         upper = ifelse(grepl("^Z", parameter), exp(-upper), upper),
         lower = ifelse(grepl("^Z", parameter), exp(-lower), lower)) %>%
ggplot(aes(year, value, color = parameter)) + geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_point() +
  facet_wrap(~parameter, scales = "free_y") + ylab("Yearly survival")


# assumptions are uncertainty of data (0.01, )

summary(res)
res2 <- do.call(optim, obj)
res2$value
res2$par




#connereturns <- filter(connereturnsdat, year %in% 1989:2016)
# saveRDS(connereturns, file = "connereturns.rds")
# connereturns <- readRDS("connereturns.rds")

# creating data list to use in rstan
# connedata <- list(#years = connereturns$year[-1],
#                   smolts = connereturns$smolt[-length(connereturns$smolt)],
#                   loggrilse = connereturns$small[-1],
#                   logSW2 = connereturns$large[-1],
#                   #N =  length(connereturns$small[-1]),
#                   cv = 0.01)