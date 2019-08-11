library(tidyverse)
library(TMB)

#compile("analysis/lc_model_hierarchical.cpp")
compile("analysis/lc_model_hierarchical.cpp", "-O0 -g") # verbose output for gdbsource

# dynamic link to object
dyn.load(dynlib("analysis/lc_model_hierarchical"))

alldata <- readRDS(file = "data/alldata.rds")

# Saving filtered data frame in new object (connereturns)
npops <- dim(alldata[[1]])[1]
nmaxyears <-  dim(alldata[[1]])[2]

# matrix(rep(-log(0.08), length(alldata$logsmolts)), nrow = npops)

pars <- list( # These first parameters are matrices
             Z1 = matrix(rep(-log(0.08), length(alldata$logsmolts)), nrow = npops), 
             Z2 = matrix(rep(-log(0.1), length(alldata$logsmolts)), nrow = npops),
             logitPr =  matrix(rep(0.1, length(alldata$logsmolts)), nrow = npops),
             logsmolts_true =  matrix(rep(mean(alldata$logsmolts, na.rm = TRUE), length(alldata$logsmolts)), nrow = npops),
             # These are population-level vectors
             logsmolts_mean = apply(alldata$logsmolts, 1, mean, na.rm = TRUE),
             logsmolts_logsd = apply(alldata$logsmolts, 1, sd, na.rm = TRUE),
             Z1_mean = rep(2.5, npops),
             Z1_logsd = rep(1.1, npops),
             Z2_mean = rep(2.5, npops),
             Z2_logsd = rep(1.1, npops),
             logitPr_mean = rep(0.1, npops),
             logitPr_logsd = rep(1.1, npops),
             # These are single "global" values
             global_Z1_mean = 2.5,
             global_Z1_logsd = 1.1,
             global_Z2_mean = 2.5,
             global_Z2_logsd = 1.1,
             global_logitPr_mean = 0.1,
             global_logitPr_logsd = 1.1
)

lapply(alldata, dim)
#mu = rep(mean(salmondata$loggrilse), length(salmondata$smolts)),
#mu2 =  rep(mean(salmondata$logSW2), length(salmondata$smolts)))

alldata["river_name"]

#alldata2 <- alldata %>% purrr::discard(names(.) %in% "river_name")
#alldata2 <- rlist::list.remove(alldata, c("river_name"))


alldata2 <- alldata %>% purrr::keep(names(.) %in% c("logsmolts", "logsmolts_SE",
                                                    "loggrilse", "logSW2",
                                                    "cv", "nyears"))
# changing cv and nyears from matrices to vectors
alldata2$cv <- as.vector(alldata2$cv)
alldata2$nyears <- as.vector(alldata2$nyears)

# removing NAs to troubleshoot model having infinite gradient
# alldata2$loggrilse[is.na(alldata2$loggrilse)] <- 6
# alldata2$logsmolts[is.na(alldata2$logsmolts)] <- 10
# alldata2$logsmolts_SE[is.na(alldata2$logsmolts_SE)] <- 2
# alldata2$logSW2[is.na(alldata2$logSW2)] <- 5

# checking the lengths of parameters, they should match expectations from TMB model
lapply(pars, dim)
lapply(alldata2, dim)
lapply(pars, length)
lapply(alldata2, length)
alldata2$nyears
alldata$river_name

# gdbsource("analysis/lc_model_hierarchical.R")

# # The model runs without any random effects or with only
# # the "top-level" random effects:
# objnore <- MakeADFun(alldata2, pars, DLL = "lc_model_hierarchical",
#                      random = c("Z1_mean", "Z2_mean", "logitPr_mean"))
# resnore <- nlminb(objnore$par, objnore$fn, objnore$gr)
# resultsnore <- sdreport(objnore)


# running hierarchical model
obj <- MakeADFun(alldata2, pars, DLL = "lc_model_hierarchical",
               random = c("logsmolts_true", "Z1", "Z2", "logitPr", 
                          "Z1_mean", "Z2_mean", "logitPr_mean"),
               control = list(eval.max = 10000, iter.max = 10000))

# This will print every parameter passed to obj$fn
obj$env$tracepar <- TRUE

str(obj)
# optimizing model objective
res <- nlminb(obj$par, obj$fn, obj$gr)
res$objective
res$par

# Calculating SDs using Laplace approximation
results <- sdreport(obj)
