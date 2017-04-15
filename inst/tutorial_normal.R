# remove all objects from R environment
rm(list = ls())
# load package
library(PET)
# fix the random seed
set.seed(17)

N <- 100000
p <- 5
mean <- rnorm(p)
sub <- t(t(rexp(p)))

library(covreg)
S <- diag(1, p, p)
for (i in 1:p){
  for (j in 1:p){
    S[i,j] <- 0.9^(abs(i-j))
  }
}
covariance <- covreg::rwish(S, 5)
# sample from multivariate Normal
x <- rmvnorm(N, mean, covariance)
summary(as.numeric(abs((cov(x) - covariance) / covariance)))

x_mvtnorm <- mvtnorm::rmvnorm(N, mean, covariance)
summary(as.numeric(abs((cov(x_mvtnorm) - covariance) / covariance)))

summary(abs((dmvnorm(x, mean, covariance) - mvtnorm::dmvnorm(x, mean, covariance, log = TRUE))/mvtnorm::dmvnorm(x, mean, covariance, log = TRUE)))

# based on Cholesky
C <- chol(covariance)
x_chol <- rmvnorm_cholesky(N, mean, C)
summary(as.numeric(abs((cov(x_chol) - covariance) / covariance)))

# based on Cholesky factor of the precision matrix (this time, apply a transpose)
C_inverse <- t(chol(solve(covariance)))
summary(abs((dmvnorm_cholesky_inverse(x, mean, C_inverse) - mvtnorm::dmvnorm(x, mean, covariance, log = TRUE))/mvtnorm::dmvnorm(x, mean, covariance, log = TRUE)))

