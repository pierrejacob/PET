# remove all objects from R environment
rm(list = ls())
# load package
library(PET)
# fix the random seed
set.seed(17)


rproposal <- function(n) rnorm(n, mean = 0, sd = 1)
dproposal <- function(x) dnorm(x, mean = 0, sd = 1, log = TRUE)

dtarget <- function(x) dnorm(x, mean = 0, sd = 0.1, log = TRUE)

N <- 10000
X <- rproposal(N)
f <- function(gamma){
  logw <- gamma * (dtarget(X) - dproposal(X))
  nw <- normalize_weight(logw)$nw
  return(1/sum(nw^2))
}
# look for the temperature ensuring an ESS of N / 2
search_gamma(0, f, 0.5*N)
