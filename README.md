# Football Task scripts
## Computational Model

Likelihood function for the computational model of the football task (Eldar et al., 2018)

USAGE: [lik, latents] = lik_football(P,data)

INPUTS:
 - P - structure of S parameter samples, with the following fields:
   - .ba1 - [S x 1] player type 1's scoring
   - .ba2 - [S x 1] player type 2's scoring
   - .ba3 - [S x 1] player type 3's scoring
   - .ba4 - [S x 1] player type 4's scoring
   - .valtype1 - [S x 1] mixes rounded and unrounded scoring for player type 1
   - .valtype2 - [S x 1] mixes rounded and unrounded scoring for player type 2
   - .valtype3 - [S x 1] mixes rounded and unrounded scoring for player type 3
   - .valtype4 - [S x 1] mixes rounded and unrounded scoring for player type 4
   - .varnocount - [S x 1] variance per minimially processed player 
   - .timehitnocountvar - [S x 1] effect of trial time on 'varnocount'
   - .varcounthit - [S x 1] decreases variance for optimally processed players
   - .w12 - [S x 1] priority weight for second ranked player type
   - .w13 - [S x 1] priority weight for third ranked player type
   - .w14 - [S x 1] priority weight for fourth ranked player type
   - .w4 - [S x 1] player type 4's priority;
   - .sequence - [S x 3] mixing coefficients for type-based, numerosity-based, and screen-location-based prioritization
   - .thresh - [S x 1] resources required for optimal processing
   - .timehitthresh - [S x 1] effect of trial time on '.thresh'
   - .inattmean - [S x 1] default scoring coefficient 
   - .counthitcapall - [S x 1] number of players beyond which default scoring is assumed
 - data - struture of experimental data with the following fields:
   - .goals - [1 x T] the participant's answer (between 0 and 10 goals)
   - .stim - [1 x 4 x T] number of players of each type
   - .dectime - [1 x T] time avalailable for deliberation (1 or 2 seconds)
   - .distance - [1 x 4 x T] average distacne from center of screen for each player type

OUTPUTS:
 - lik - [S x 1] log-likelihoods
 - latents - a structure with the following fields:
   - .goals - [S x T] model's random answer
   - .goals_max - [S x T] model's most likely answer
   - .processing_weights - [S x 4 x T] model's resource allocation
## Sequenceness analysis
Compute sequenceness between pairs of time series

USAGE: seq = sequenceness(X, wind, maxgap)

INPUTS: 
 - X - [T x Q x S] data matrix containing S time series for T trials and Q timepoints per trial
 - wind - length of time window to use for each calculation of sequenceness 
 - maxgap - maximal time lag between time series to consider

OUTPUTS:
 - seq - [maxgap x T x P x Q-wind] sequencesness for each time lag upto maxgap, for each trial, for each pair of time series, for each starting timepoints
