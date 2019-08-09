Hierarchical life cycle model to estimate marine survival in Atlantic salmon
=====

This repository contains the source code for the TMB model estimating marine survival of
Atlantic salmon using time series abundance data of returning adults and outmigrating
smolts. 

The model is defined in `analysis/lc_model_hierarchical.cpp`.

The non-hierarchical version of the model (i.e. where survival is estimates for individual
populations) is defined in  `analysis/lc_model.cpp`.

So far, the rivers for we have sufficient data to run the model for are :

- Nashwaak River, NB (St John River)
- Conne River, NL 
- Saint-Jean River, QC
- Trinité River, QC

These data are loaded from separate repositories but merged into .rds files for easier use.

At this point the non-hierarchical model running separately for each river works, but the 
hierarchical one is not converging w.
I'm also hoping to incorporate a stock-recruitment function to further inform *smolts_true* and 
also to estimate egg-to-smolt survival.



