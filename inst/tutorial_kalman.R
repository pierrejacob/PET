rm(list = ls())
set.seed(17)

## Model
## X_0 ~ Normal(mu_0, Sigma_0)
## X_t = Phi X_{t-1} + Sigma_W W_t
## Y_t = A X_t + Sigma_V V_t

Phi <- matrix(c(1,1,0,0,
                0,1,0,0,
                0,0,1,1,
                0,0,0,1), byrow = TRUE, ncol = 4)
A <- matrix(c(1,0,0,0,
              0,0,1,0), byrow = TRUE, ncol = 4)

T <- 100
X <- matrix(0, ncol = 4, nrow = T+1)
Y <- matrix(0, ncol = 2, nrow = T)
sigma_W <- 0.1
sigma_V <- 1
for (t in 1:100){
  X[t+1,] <- Phi %*% matrix(X[t,], ncol = 1) + rnorm(4, 0, sigma_W)
  Y[t,] <- A %*% matrix(X[t+1,], ncol = 1) + rnorm(2, 0, sigma_V)
}

plot(Y[,1], Y[,2])
lines(X[,1], X[,3])

### Kalman filter
kf_means <- matrix(0, ncol = 4, nrow = T+1)
m_current <- rep(0, 4)
kf_means[1,] <- m_current
V_current <- matrix(0, ncol = 4, nrow = 4)
Sigma_W <- diag(sigma_W^2, 4, 4)
Sigma_V <- diag(sigma_V^2, 2, 2)
Id <- diag(1, 4, 4)
for (t in 1:100){
  # prediction step
  m_next <- Phi %*% matrix(m_current, ncol = 1)
  V_next <- Phi %*% V_current %*% t(Phi) + Sigma_W
  # update step
  K <- V_next %*% t(A) %*% solve(A %*% V_next %*% t(A) + Sigma_V)
  m_current <- m_next + K %*% (Y[t,] - A %*% m_next)
  V_current <- (Id - K %*% A) %*% V_next
  kf_means[t+1,] <- m_current
}

lines(kf_means[,1], kf_means[,3], col = "blue")

# compare to existing package
library(astsa)
kfastsa <- Kfilter0(T, Y, A, rep(0, 4), diag(0, 4, 4), Phi, chol(Sigma_W), chol(Sigma_V))
lines(kfastsa$xf[1,1,], kfastsa$xf[3,1,], col = "orange", lty = 3)

### Kalman smoother
# requires storing more stuff in the forward pass
kf_means <- matrix(0, ncol = 4, nrow = T+1)
kf_predmeans <- matrix(0, ncol = 4, nrow = T)

kf_variances <- list()
kf_predvariances <- list()
m_current <- rep(0, 4)
kf_means[1,] <- m_current
V_current <- matrix(0, ncol = 4, nrow = 4)
kf_variances[[1]] <- V_current
Sigma_W <- diag(sigma_W^2, 4, 4)
Sigma_V <- diag(sigma_V^2, 2, 2)
Id <- diag(1, 4, 4)
for (t in 1:100){
  # prediction step
  m_next <- Phi %*% matrix(m_current, ncol = 1)
  V_next <- Phi %*% V_current %*% t(Phi) + Sigma_W
  kf_predmeans[t,] <- m_next
  kf_predvariances[[t]] <- V_next
  # update step
  K <- V_next %*% t(A) %*% solve(A %*% V_next %*% t(A) + Sigma_V)
  m_current <- m_next + K %*% (Y[t,] - A %*% m_next)
  V_current <- (Id - K %*% A) %*% V_next
  kf_means[t+1,] <- m_current
  kf_variances[[t+1]] <- V_current
}


# smoothing backward pass
ks_means <- matrix(0, ncol = 4, nrow = T+1)
ks_means[T+1,] <- kf_means[T+1,]
ks_variances <- list()
ks_variances[[T+1]] <- kf_variances[[T+1]]
for (t in 99:1){
  L <- kf_variances[[t+1]] %*% t(Phi) %*% solve(kf_predvariances[[t+1]])
  m <- kf_means[t+1,] + L %*% (ks_means[t+2,] - kf_predmeans[t+1,])
  V <- kf_variances[[t+1]] + L %*% (ks_variances[[t+2]] - kf_predvariances[[t+1]])
  ks_means[t+1,] <- m
  ks_variances[[t+1]] <- V
}

lines(ks_means[,1], ks_means[,3], col = "red")
library(astsa)
ksastsa <- Ksmooth0(T, Y, A, rep(0, 4), diag(0, 4, 4), Phi, chol(Sigma_W), chol(Sigma_V))
lines(ksastsa$xs[1,1,], ksastsa$xs[3,1,], col = "orange", lty = 3)

