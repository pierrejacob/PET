# remove all objects from R environment
rm(list = ls())
# load package
library(PET)
# fix the random seed
set.seed(17)

N <- 1000
logweights <- rnorm(N)

normalize_weight_results <- normalize_weight(logweights)
normalized_weights <- normalize_weight_results$nw

# Now resampling schemes, coded in R
N <- 6
logweights <- rnorm(N)
normalize_weight_results <- normalize_weight(logweights)
normalized_weights <- normalize_weight_results$nw

Nprime <- 10000
ancestors <- multinomial_resampling_R(Nprime, normalized_weights)
summary(abs((tabulate(ancestors) / Nprime - normalized_weights)/normalized_weights))

ancestors <- multinomial_resampling(Nprime, normalized_weights)
summary(abs((tabulate(ancestors) / Nprime - normalized_weights)/normalized_weights))

ancestors <- systematic_resampling(Nprime, normalized_weights)
summary(abs((tabulate(ancestors) / Nprime - normalized_weights)/normalized_weights))

ancestors <- replicate(Nprime,tabulate(ssp_resampling(normalized_weights)))
summary(abs((rowMeans(ancestors/N)-normalized_weights)/normalized_weights))

