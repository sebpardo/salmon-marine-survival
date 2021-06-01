//  I would like to suggest that you do some quick sensitivity tests: 
// (1) use a relatively weak informative prior for Z (equations 12 and 13), e.g., using sigma = 0.5 or 0.8; 
// (2) increase/decrease the median values for Z by 10% or 20%; 
// (3) you may also do the similar tests for parameter Pg,r,t in equations 14-16. 

// RStan model of marine survival that is simplified, with the only hierarchical part being 
// the estimation of Pr which is not logit transformed but probit transformed
// Also, the observation error in smolts and hierarchical estimation of Pr are reparameterized
// to be non-centered in order to avoid any funneling/sampling space problems.

data {
  int<lower=1> N;         // number of populations
  int<lower=1> K;         // total number of observations
  real z1_mu;          // Mu prior for Z1 
  real z2_mu;          // Mu prior for Z2
  real<lower=0> z1_sd; // SD prior for Z1 
  real<lower=0> z2_sd; // SD prior for Z2

  // instead of two-dimensional arrays, in Stan these need to be coded as single long vectors
  real logsmolts[K];      //
  real<lower=0> logsmolts_cv[K];   //
  real loggrilse[K];      //
  real logSW2[K];         //
  
  real logSW1_cv[K];       // annual estimates of CV of log transformed 1SW returns
  real logSW2_cv[K];       // annual estimates of CV of log transformed 2SW returns

  // one-dimensional
 // vector<lower=0>[N] returns_cv;   // coefficient of variation of 1SW and 2SW returns
  int nyears[N];          // number of years of data for each population, used for indexing
  real logis_mu[N];       // mu parameter of logistic distribution for Pr prior for each population
  real logis_sigma[N];    // sigma parameter of logistic distribution for Pr prior for each population
}

// this are the parameters for which we do sampling (i.e. assign priors)
parameters {
  // two-dimensional arrays
  real<lower=0> Z1[K];   //
  real<lower=0> Z2[K];   //
  vector[K] beta_raw;   // scaling factor of the hierarchical logitPr estimates
  vector[K] logsmolts_raw;   // scaling factor of the observation error hierarchical estimates
  // population level
  real logsmolts_mean[N];
  vector<lower=0>[N] logsmolts_sigma;
  real logitPr_mu[N];
  vector<lower=0>[N] logitPr_sigma;
}

// estimate additional parameters from sampled parameters (don't have priors)
transformed parameters {  
  vector[K] logitPr;
  vector[K] logsmolts_true;   //
  real<lower=0, upper=1> Pr[K];   //
  real<lower=0, upper=1> Pr_mu[N];   //

  { 
    //starting BLOCK, objects in here cannot be seen outside (this is a
    // trick for declaring integers in the transformed parameters section)
    int pos = 1;
  
    for (j in 1:N) {
        int last = pos + nyears[j] - 1;
        // indexing in left side of assignment is just like R's
        logitPr[pos:last] =  logitPr_mu[j] + segment(beta_raw, pos, nyears[j]) * logitPr_sigma[j]; 
        logsmolts_true[pos:last] = logsmolts_mean[j] + segment(logsmolts_raw, pos, nyears[j]) * logsmolts_sigma[j];

        Pr_mu[j] = Phi_approx(logitPr_mu[j]); // either Phi, Phi_approx, of inv_logit
        // inverse logit link function for Pr
        // Pr_mu[j] = exp(logitPr_mu[j])/(1 + exp(logitPr_mu[j]));
        pos = pos + nyears[j]; 
    }
  }
  
  for (i in 1:K) {
    // Pr[i] = exp(logitPr[i])/(1 + exp(logitPr[i])); // inverse logit link function for Pr
    Pr[i] = Phi_approx(logitPr[i]);
  }
}

// parameters declared here are not saved to output
model {
  int pos2;      // indexing value that is updated for extracting different segment() 
  real mu[K];   // estimated log 1SW 
  real mu2[K];  // estimated log 2SW
  pos2 = 1;

// these estimates of mu and mu2 are additive a they are in log
// space, thus proportion survival S1 has been transformed to 
// intstantaneous total mortality Z1
  for (i in 1:K) {
    mu[i] = logsmolts_true[i] + log(Pr[i]) - Z1[i]; // estimated log1SW
    mu2[i] = logsmolts_true[i] - Z1[i] + log1m(Pr[i])  - Z2[i]; // estimated log2sw
  }
  
// Because Stan can't deal with ragged arrays/matrices and we have a different number of years
// for each population, the matrices are flattened into a single 
// vector and the indexing is done by separately looping through each "chunk" of the vector
// using the segment() function

  for (j in 1:N) {
      // observation error in smolts (dist of error around x_true)
      segment(logsmolts, pos2, nyears[j]) ~ normal(segment(logsmolts_true, pos2, nyears[j]), logsmolts_cv[j]); 
 
      // non-centered parameterizations
      segment(beta_raw, pos2, nyears[j]) ~ std_normal();  // implies logitPr ~ normal(logitPr_mu, logitPr_sigma)
      segment(logsmolts_raw, pos2, nyears[j]) ~ std_normal(); // implies logsmolts_true ~ normal(logsmolts_mean, logsmolts_sd)
      // obs error in 1SW and 2SW (based on data cv) log transformed SD is roughly equal to log transformed CV
      
      segment(loggrilse, pos2, nyears[j]) ~ normal(segment(mu, pos2, nyears[j]), segment(logSW1_cv, pos2, nyears[j])); 
      segment(logSW2, pos2, nyears[j]) ~ normal(segment(mu2, pos2, nyears[j]), segment(logSW2_cv, pos2, nyears[j]));
      
      pos2 = pos2 + nyears[j];
    }
    
    // Priors
    logitPr_mu ~ normal(logis_mu, logis_sigma);   // prior of Pr for each population
    logitPr_sigma ~ std_normal();                 // uninformative prior for logitPr_sigma
    Z1 ~ lognormal(z1_mu, z1_sd);    
    Z2 ~ lognormal(z2_mu, z2_sd);
}

generated quantities {
  real<lower=0, upper=1> S1[K];
  real<lower=0, upper=1> S2[K];
  real<lower=0> smolts_true[K];
  for (i in 1:K) {
    S1[i] = exp(-Z1[i]);
    S2[i] = exp(-Z2[i]);
    smolts_true[i] = exp(logsmolts_true[i]);
  }
}
