// Simple linear regression.
#include <TMB.hpp>
template<class Type>
// TMB uses an objective function that returns the negative log-likelihood
// which is minimized in R using an optimization function (e.g. nlminb)
Type objective_function<Type>::operator() ()
{
  // Data section (needs to be log transformed before being loaded)
  DATA_ARRAY(logsmolts); // yearly smolt abundance in log space
  DATA_ARRAY(logsmolts_SE); // standard error of smolt abundance
  DATA_ARRAY(loggrilse); // yearly grilse abundance in log space
  DATA_ARRAY(logSW2);    // yearly 2SW abundance in log space
  //DATA_VECTOR(Pr); // Simplifying the model by including Pr as data instead of parameter
  DATA_VECTOR(cv);        // CV of adult returns estimates (observation error) for each population
  DATA_VECTOR(nyears); // number of years in each population, needed for indexing
  
  // int npop = logsmolts.size(); // number of populations (nrow)
  // int n = logsmolts[0].size(); // max number of years (ncol) not needed for indexing
  //int npop = logsmolts.size()/cv.size();
  //int n = cv.size();
  
  int npop = logsmolts.rows();
  int n = logsmolts.cols();
  
  // Yearly parameters across all populations
  PARAMETER_ARRAY(Z1); // yearly instantaneous mortality rate of first year at sea
  PARAMETER_ARRAY(Z2); // yearly instantaneous mortality rate of additional (second) year at sea
  PARAMETER_ARRAY(logitPr); // yearly proportion returning as grilse in logit 
  PARAMETER_ARRAY(logsmolts_true); // True yearly smolt abundance in log space
  // distribution parameters (vector, as one for each population)
  PARAMETER_VECTOR(logsmolts_mean); // mean log smolt abundance across years
  PARAMETER_VECTOR(logsmolts_logsd); // sd log smolt abundance across years
  PARAMETER_VECTOR(Z1_mean); // mean global mortality 
  PARAMETER_VECTOR(Z1_logsd); // global mortality st dev
  PARAMETER_VECTOR(Z2_mean); // mean 
  PARAMETER_VECTOR(Z2_logsd); //
  PARAMETER_VECTOR(logitPr_mean); // 
  PARAMETER_VECTOR(logitPr_logsd); //
  // global parameters (single values)
  PARAMETER(global_Z1_mean); // mean global mortality 
  PARAMETER(global_Z1_logsd); // global mortality st dev
  PARAMETER(global_Z2_mean); // mean 
  PARAMETER(global_Z2_logsd); //
  PARAMETER(global_logitPr_mean); // 
  PARAMETER(global_logitPr_logsd); //
  // 
  //PARAMETER_VECTOR(mu); // estimated returns of grilse
  //PARAMETER_VECTOR(mu2); // estimated 2SW returns
  // Log Likelihood
  Type nll=0; 
  //Matrix<typename double, int npop, int n> 
  array<Type> Pr(npop, n); // need to double check that this indexing is correct
  array<Type> mu(npop, n); // seems that the format is indeed (columns, rows)
  array<Type> mu2(npop, n);
  
  // first index of elements in C++ is 0, not 1 like R
  for (int j = 0; j < npop; j++) {  // first loop is for populations
  
    for (int i = 0; i < nyears(j); i++) { // second loop is for years, different for each popn
      
      Pr(j,i) = exp(logitPr(j,i))/(1 + exp(logitPr(j,i))); // logit link 
  
      // these estimates of mu and mu2 are additive a they are in log
      // space, thus proportion survival S1 has been transformed to 
      // intstantaneous total mortality Z1
      mu(j,i) = logsmolts(j,i) + log(Pr(j,i)) - Z1(j,i); // estimated loggrilse
      mu2(j,i) = logsmolts(j,i) + log(1 - Pr(j,i)) - Z1(j,i) - Z2(j,i); // estimated log2sw
      
      nll -= dnorm(Z1(j,i), Z1_mean(j), exp(Z1_logsd(j)), true);
      nll -= dnorm(Z2(j,i), Z2_mean(j), exp(Z2_logsd(j)), true);
      
      nll -= dnorm(logitPr(j,i), logitPr_mean(j), exp(logitPr_logsd(j)), true);
      
      nll -= dnorm(logsmolts_true(j,i), logsmolts_mean(j), exp(logsmolts_logsd(j)), true); // distribution of x_true
      nll -= dnorm(logsmolts(j,i), logsmolts_true(j,i), logsmolts_SE(j,i), true); // oibservation error in x (dist of error around x_true)
      
      nll -= dnorm(loggrilse(j,i), mu(j,i), cv(j), true); // observation error in grilse (based on data cv)
      nll -= dnorm(logSW2(j,i), mu2(j,i), cv(j), true); // observation error in 2SW (based on data cv)
      
      // NLL DECLARATION FOR GLOBAL PARAMETER ESTIMATES
      nll -= dnorm(Z1_mean(j), global_Z1_mean, exp(global_Z1_logsd), true);
      nll -= dnorm(Z2_mean(j), global_Z2_mean, exp(global_Z2_logsd), true);
      nll -= dnorm(logitPr_mean(j), global_logitPr_mean, exp(global_logitPr_logsd), true);
  }
}
  ADREPORT(Pr);

  return(nll);
}
