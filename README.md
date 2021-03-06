Hierarchical life cycle model to estimate marine survival in Atlantic salmon
=====

This repository contains the source code for the model estimating marine
survival of Atlantic salmon using time series abundance data of returning
adults and outmigrating smolts, written both in Stan and TMB (although the Stan
version is the most up-to-date and the one used for publication). 

The Stan model used in for publication defined in `analysis/survival_hier_ncp_annualcvs.stan`. 
The TMB model is defined in
`analysis/lc_model_hierarchical.cpp` but was not actively developed; it is left
here as a reference.

The non-hierarchical version of the model (i.e. where survival is estimates for
individual populations) is defined in  `analysis/lc_model.cpp`.

So far, the rivers for we have sufficient data to run the model for are:

- LaHave River, NS
- Nashwaak River, NB (St John River)
- Conne River, NL 
- Campbellton River, NL 
- Western Arm Brook, NL 
- Saint-Jean River, QC
- Trinité River, QC

These data are loaded from separate repositories but merged into .rds files for
easier use.

The model is hierarchical in the estimation of *Pr*, where the yearly value of
this parameter is estimated based on a population-level distribution, which has
relatively informative priors for each river based on their known life
histories, and also in the estimation of true smolt abundance, which is done
ona population-level basis.

<!-- I'm also hoping to incorporate a stock-recruitment function to further
inform *smolts_true* and also to estimate egg-to-smolt survival.-->



