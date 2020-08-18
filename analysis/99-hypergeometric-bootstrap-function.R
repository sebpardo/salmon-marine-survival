### Function to estimate coefficient of variation for 1SW and 2SW based on 
### the total number of estimated small and large returns as well as a
### subset of the scale data to assign small and large as  1SW and 2SW returns.
### This function uses a hypergeometric distribution 

hyper_boot <- function(year = NA, n_small, n_large, scales_small_1SW, scales_small_2SW,
                       scales_small_other,
                       scales_large_1SW, scales_large_2SW, 
                       scales_large_other,
                       iter = 10000, 
                       debug = FALSE,
                       list = FALSE,
                       ...) {
  
  #if (!is.na(year) & year == 1995) browser()
  
  scales_small_total <- scales_small_1SW + scales_small_2SW + scales_small_other # total small scales samples
  scales_large_total <- scales_large_1SW + scales_large_2SW + scales_large_other  # total large scales samples
  
  if (scales_small_total > n_small) stop(paste("There are more small scale samples than returns for year", year))
  if (scales_large_total > n_large) stop(paste("There are more large scale samples than returns for year", year))
  
  # vectors of all potential abundances for small and large returns
  if (n_small > 0) {vec_nsmall <- 0:n_small} else {vec_nsmall <- 0}
  if (n_large > 0) {vec_nlarge <- 0:n_large} else {vec_nlarge <- 0}
  
  
  prob_sm1 <- dhyper(scales_small_1SW, vec_nsmall, n_small - vec_nsmall, scales_small_total) # prob of small 1SW
  prob_sm2 <- dhyper(scales_small_2SW, vec_nsmall, n_small - vec_nsmall, scales_small_total) # prob of small 2SW
  prob_la1 <- dhyper(scales_large_1SW, vec_nlarge, n_large - vec_nlarge, scales_large_total) # prob of large 1SW
  prob_la2 <- dhyper(scales_large_2SW, vec_nlarge, n_large - vec_nlarge, scales_large_total) # prob of large 2SW

  prob_la1[is.nan(prob_la1)] <- 1
  prob_la2[is.nan(prob_la2)] <- 1
  
  # Sampling vectors of abundances based on probabilities of hypergeometric distribution
  hyp_draws_sm1 <- base::sample(vec_nsmall, iter, replace = TRUE, prob = prob_sm1)
  hyp_draws_sm2 <- base::sample(vec_nsmall, iter, replace = TRUE, prob = prob_sm2)
  hyp_draws_la1 <- base::sample(vec_nlarge, iter, replace = TRUE, prob = prob_la1)
  hyp_draws_la2 <- base::sample(vec_nlarge, iter, replace = TRUE, prob = prob_la2)
  
  # adding draws of small and large 1SW and 2SW abundances, respectively, and log-transforming
  log_1SW <- log(hyp_draws_sm1 + hyp_draws_la1) # adding 1SW abundances of small + large
  log_2SW <- log(hyp_draws_sm2 + hyp_draws_la2) # adding 2SW abundances of small + large
  
  if (list) {
    return(list(year = year,
                meanlog1SW =  mean(log_1SW),
                sdlog1SW = sd(log_1SW),
                meanlog2SW = mean(log_2SW), 
                sdlog2SW = sd(log_2SW),
                sm1 = hyp_draws_sm1,
                sm2 = hyp_draws_sm2,
                la1 = hyp_draws_la1,
                la2 = hyp_draws_la2,
                prob_sm1 = prob_sm1,
                prob_sm2 = prob_sm2,
                prob_la1 = prob_la1,
                prob_la2 = prob_la2))
  } else {
    return(data.frame(year = year,
                      meanlog1SW =  mean(log_1SW),
                      sdlog1SW = sd(log_1SW),
                      meanlog2SW = mean(log_2SW), 
                      sdlog2SW = sd(log_2SW)))    
  }
}


hyper_boot(n_small = 511,
           n_large = 131,
           scales_small_1SW = 500, 
           scales_small_2SW = 10,
           scales_small_other = 1,
           scales_large_1SW = 30, 
           scales_large_2SW = 100, 
           scales_large_other = 1,
           #debug = TRUE,
           iter = 10000)


