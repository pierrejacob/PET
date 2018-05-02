# remove all objects from R environment
rm(list = ls())
# load package
library(PET)
# fix the random seed
set.seed(17)
# load module
module_tree <<- Module("module_tree", PACKAGE = "PET")
# TreeClass is a class of trees
TreeClass <<- module_tree$Tree

# we illustrate the use of the class via a simple particle filter algorithm storing the surviving trajectories
N <- 100 # number of particles
T <- 10 # number of observations
observations <- rnorm(T)
# model:
# x_1 ~ Normal(0,1); x_t | x_{t-1} ~ Normal(rho x_{t-1}, 1) where rho = 0.9
# y_t | x_t ~ Normal(x_t, 1)
# create new instance of the class
Tree <- new(TreeClass, N, 10*N, 1) # first argument is number of leaves; second is expected total size of tree; third is dimension of each leaf
# first particles
x <- rnorm(N, mean = 0, sd = 1)
# initialize tree
Tree$init(matrix(x, nrow = 1, ncol = N))
weights <- dnorm(x = observations[1], mean = x, sd = 1)
weights <- weights / sum(weights)
ancestors <- sample(x = 1:N, size = N, replace = T, prob = weights)
for (t in 2:T){
  x <- 0.9 * x[ancestors] + rnorm(N, mean = 0, sd = 1)
  weights <- dnorm(x = observations[1], mean = x, sd = 1)
  weights <- weights / sum(weights)
  Tree$update(matrix(x, nrow = 1, ncol = N), ancestors - 1) # note that Tree functions are in C++ and so follow the C++ conventions on indices
  ancestors <- sample(x = 1:N, size = N, replace = T, prob = weights)
}
# Now we can query the tree, for instance, getting all the trajectories
paths <- matrix(0, nrow = N, ncol = T)
for (i in 1:N){
  paths[i,] <- Tree$get_path(i-1)
}
# plot all the trajectories
matplot(t(paths), type = "l", col = rgb(0,0,0, alpha = 0.2), lty = 1)
title("Paths (in C++)")

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# Same task as previously using a tree implemented in R with no dependence on
# the C++ Tree module

# fix the random seed
set.seed(17)
# we illustrate the use of the class via a simple particle filter algorithm storing the surviving trajectories
N <- 100 # number of particles
T <- 10 # number of observations
observations <- rnorm(T)
# model:
# x_1 ~ Normal(0,1); x_t | x_{t-1} ~ Normal(rho x_{t-1}, 1) where rho = 0.9
# y_t | x_t ~ Normal(x_t, 1)
# first particles
# WARNING: particles must be formated as dimX by N matrices
x <- matrix(rnorm(N, mean = 0, sd = 1), ncol = N)
# initialize the tree
Tree_R = tree_init(x)
# reweight and resample
weights <- dnorm(x = observations[1], mean = x, sd = 1)
weights <- weights / sum(weights)
ancestors <- sample(x = 1:N, size = N, replace = T, prob = weights)
for (t in 2:T){
  x <- matrix(0.9 * x[ancestors] + rnorm(N, mean = 0, sd = 1), ncol = N)
  weights <- dnorm(x = observations[1], mean = x, sd = 1)
  weights <- weights / sum(weights)
  Tree_R = tree_update(Tree_R, x, ancestors)
  ancestors <- sample(x = 1:N, size = N, replace = T, prob = weights)
}
# Now we can query the tree, for instance, getting all the trajectories
paths_R <- matrix(0, nrow = N, ncol = T)
for (i in 1:N){
  paths_R[i,] <- tree_getpath(Tree_R,i)
}
# plot all the trajectories
matplot(t(paths_R), type = "l", col = rgb(0,0,0, alpha = 0.2), lty = 1)
title("Paths (in R)")

# Check that both implementations give the same paths
all(paths == paths_R)
