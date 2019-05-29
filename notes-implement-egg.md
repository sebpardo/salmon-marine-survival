February 13 2019, discussion with Geir
========

## Accounting for differential exposure in S1 for 1SW and 2SW maidens

One big assumption of estimating S1 is that the value of S1 is the same for both returning grilse and non-returning 2SW. However, we know that is not true given that 
grilse leave wintering grounds earlier and spend some months travelling back to their natal rivers. We can account for this by:

- We can multiply S1 by proportion of number of months that grilse spend at sea vs that of those who dont migrate to spawn
- how much of the mortality happens on the back migration (proportion)

## Implementing egg production into model

The idea is to feedback data from returns into the model to inform smolt estimates based on egg production (and also estimate egg-to-smolt survival):

returns -> spawners -> eggs -> smolts

This would result in a "true" life cycle model, and would 1) improve the estimation of smolts_true by "grounding" expected values based on eggs, and 2) produce estimates of egg-to-smolt survival resulting from partial pooling (by sampling yearly values from a distribution).

Assumptions:

- error is normally distributed, accounts for annual variation in eggs to smolt survival
- river age is constant (3) thus eggs all come from a single cohort (although they could be lagged among spawner cohorts but we shouldn't do this to begin with)

instead of 

$$logsmolts_{true,t} \sim N(\mu_{smolts}, \sigma_{smolts})$$

we do this 

$$logsmolts_{true,t} \sim N(f(eggs_{t-4}), \sigma_{smolts})$$

with a simple function to estimate egg production:

$$f(eggs) = log(eggs_true) - a$$

or even better, 

f(eggs) = a + b * logeggs_{t-4}

with

$logeggs_{true} \sim N(f(1SW,2SW), \sigma_{eggs})$

Different ways of doing thesame thing: Estimating f(eggs):

$smolts_{true,t}= (a + b * logeggs_{t-4}) + \epsilon_t$ <!-- bracket is f(eggs) -->
$logsmolts_t \sim N(smolts_{true,t}, exp(logsmolts_{logsd}))$
$\epsilon_t \sim N(0, \sigma_smolts)$ <!-- epsilon is the variance around expected egg-to-smolt survival !!!-->

OR

$logsmolts_t \sim N(a+b+*logeggs{t-4} + \epsilon_t, exp(logsmolts_{logsd}))$
$\epsilon_t \sim N(0, \sigma_smolts)$

where mt and et are the measurement and the variance respectively. $\epsilon_t$ captures the variance in egg-to-smolt survival, which have the following distributions:

$logsmolts_{true,t} = (a + b * logeggs_{t-4}) + \epsilon_t$

```
DATA_VECTOR(EGGS)

PARAMETER_VECTORS(EPSILON)
PARAMETER_VECTORS(M)
PARAMETER(A)
PARAMETER(B)
PARAMTER(SIMGA_SMOLTS)
```

-------

$$grilse_t = S1_t * Pr_t* smolts_{t-1}$$

$$2SW_{t+1} = S1_t * K_{exposure} * (1 - Pr_t) * S2_{t+1} * smolts_{t-1}$$

------

![](geir-model-notes.jpg)
 
