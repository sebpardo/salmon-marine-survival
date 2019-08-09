# Merge data from all populations into structure to input to TMB hierachical model

library(tidyverse)
library(rlist)
#library(data.table)
library(rowr)

# NEED TO CHECK FOR -Inf VALUES!!!

saintjeandata <- readRDS(file = "data/saintjeandata.rds")
trinitedata <- readRDS(file = "data/trinitedata.rds")
connedata <- readRDS(file = "data/connedata.rds")
nashwaakdata <- readRDS(file = "data/nashwaakdata.rds")

# Converting lists for separate population into single list with  matrices for each data type 
# where each population is a separate row
# based on answer in 
# https://stackoverflow.com/questions/57401290/bind-vectors-across-lists-to-single-list-of-matrices

alldata <- list(saintjeandata, trinitedata, connedata, nashwaakdata) %>% 
  purrr::transpose() %>%
  map(~ reduce(.x, cbind.fill, fill = NA) %>% 
        t %>% 
        `row.names<-`(NULL))


saveRDS(alldata, file = "data/alldata.rds")
