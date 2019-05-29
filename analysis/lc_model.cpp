// Simple linear regression.
#include <TMB.hpp>
template<class Type>
// TMB uses an objective function that returns the negative log-likelihood
// which is minimized in R using an optimization function (e.g. nlminb)
Type objective_function<Type>::operator() ()
{
  // Data section (needs to be log transformed before being loaded)
  DATA_VECTOR(logsmolts); // yearly smolt abundance in log space
  DATA_VECTOR(logsmolts_SE); // standard error of smolt abundance
  DATA_VECTOR(loggrilse); // yearly grilse abundance in log space
  DATA_VECTOR(logSW2);    // yearly 2SW abundance in log space
  DATA_SCALAR(cv);        // CV of adult returns estimates (observation error)
  //DATA_VECTOR(Pr); // Simplifying the model by including Pr as data instead of parameter
  
  int n = logsmolts.size();
  
  // Yearly parameters
  PARAMETER_VECTOR(Z1); // yearly instantaneous mortality rate of first year at sea
  PARAMETER_VECTOR(Z2); // yearly instantaneous mortality rate of additional (second) year at sea
  PARAMETER_VECTOR(logitPr); // yearly proportion returning as grilse in logit 
  PARAMETER_VECTOR(logsmolts_true); // True yearly smolt abundance in log space
  // distribution parameters (single values)
  PARAMETER(logsmolts_mean); // mean log smolt abundance across years
  PARAMETER(logsmolts_logsd); // sd log smolt abundance across years
  PARAMETER(Z1_mean); // 
  PARAMETER(Z1_logsd); //
  PARAMETER(Z2_mean); // 
  PARAMETER(Z2_logsd); //
  PARAMETER(logitPr_mean); // 
  PARAMETER(logitPr_logsd); //
  // 
  //PARAMETER_VECTOR(mu); // estimated returns of grilse
  //PARAMETER_VECTOR(mu2); // estimated 2SW returns
  // Log Likelihood
  Type nll=0; 
  vector<Type> Pr(n);
  vector<Type> mu(n);
  vector<Type> mu2(n);
  
  // first index of elements in C++ is 0, not 1 like R
  for (int i = 0; i < n; i++) {
    
    Pr(i) = exp(logitPr(i))/(1 + exp(logitPr(i))); // logit link 

    // these estimates of mu and mu2 are additive a they are in log
    // space, thus proportion survival S1 has been transformed to 
    // intstantaneous total mortality Z1
    mu(i) = logsmolts(i) + log(Pr(i)) - Z1(i); // estimated loggrilse
    mu2(i) = logsmolts(i) + log(1 - Pr(i)) - Z1(i) - Z2(i); // estimated log2sw
    
    nll -= dnorm(Z1(i), Z1_mean, exp(Z1_logsd), true);
    nll -= dnorm(Z2(i), Z2_mean, exp(Z2_logsd), true);
    
    nll -= dnorm(logitPr(i), logitPr_mean, exp(logitPr_logsd), true);
    
    nll -= dnorm(logsmolts_true(i), logsmolts_mean, exp(logsmolts_logsd), true); // distribution of x_true
    nll -= dnorm(logsmolts(i), logsmolts_true(i), logsmolts_SE(i), true); // oibservation error in x (dist of error around x_true)
    
    nll -= dnorm(loggrilse(i), mu(i), cv, true); // observation error in grilse (based on data cv)
    nll -= dnorm(logSW2(i), mu2(i), cv, true); // observation error in 2SW (based on data cv)
  }

  ADREPORT(Pr);
  
  return(nll);
  
}
