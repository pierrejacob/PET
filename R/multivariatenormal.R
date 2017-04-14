#'@rdname rmvnorm
#'@title Generate multivariate Normals
#'@export
rmvnorm <- function(n, mean, covariance){
  return(rmvnorm_(n, mean, covariance))
}

#'@export
rmvnorm_cholesky <- function(n, mean, cholesky){
  return(rmvnorm_cholesky_(n, mean, cholesky))
}

#'@export
dmvnorm <- function(x, mean, covariance){
  return(dmvnorm_(x, mean, covariance))
}

#'@export
dmvnorm_cholesky_inverse <- function(x, mean, cholesky_inverse){
  return(dmvnorm_cholesky_inverse_(x, mean, cholesky_inverse))
}

