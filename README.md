Hierarchical life cycle model to estimate marine survival in Atlantic salmon
=====

This repository contains the source code for the TMB model estimating marine survival of
Atlantic salmon using time series abundance data of returning adults and outmigrating
smolts. 

The model is defined in `lc_model.cpp`. 

So far, the rivers for which the model is ran are:

- Nashwaak River, NB (St John River)
- Conne River, NL (**still need to transform small/large into 1SW/2SW/RS**)

These data are loaded from separate repositories.

At this point the model is ran separately for each river, but the idea is to do it hierarchically.
Also, I'm hoping to incorporate a stock-recruitment function to further inform *smolts_true* and 
also to estimate egg-to-smolt survival.



