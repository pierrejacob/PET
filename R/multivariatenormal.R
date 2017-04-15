#'@rdname rmvnorm
#'@title Generate multivariate Normals
#'@description Function to generate draws from a multivariate Normal distribution.
#'@details  This function does not check anything (i.e. that the given
#' covariance is PSD). Thus it is faster than functions in standard packages, but more risky.
#'@param n is the number of desired samples
#'@param mean is the mean vector
#'@param covariance is the covariance matrix
#'@return a n x d matrix where d is both the length of \code{mean} and the dimension of \code{covariance}. Each row contains a draw from
#'the desired multivariate Normal.
#'@examples
#' rmvnorm(10, rep(0, 5), diag(1, 5, 5))
#'@seealso \code{\link{rmvnorm_cholesky}}, \code{\link{dmvnorm}}, \code{\link{dmvnorm_cholesky_inverse}}
#'@export
rmvnorm <- function(n, mean, covariance){
  return(rmvnorm_(n, mean, covariance))
}

#'@rdname rmvnorm_cholesky
#'@title Generate multivariate Normals
#'@description Function to generate draws from a multivariate Normal distribution, given Cholesky factor of the covariance.
#'@details  This function does not check anything (i.e. that the given
#' covariance is PSD). Thus it is faster than functions in standard packages, but more risky.
#'@param n is the number of desired samples
#'@param mean is the mean vector
#'@param cholesky is the cholesky factor of the covariance matrix, obtained e.g. with \code{chol()}.
#'@return a n x d matrix where d is both the length of \code{mean} and the dimension of \code{cholesky}. Each row contains a draw from
#'the desired multivariate Normal.
#'@examples
#' rmvnorm_cholesky(10, rep(0, 5), chol(diag(1, 5, 5)))
#'@seealso \code{\link{rmvnorm}}, \code{\link{dmvnorm}}, \code{\link{dmvnorm_cholesky_inverse}}
#'@export
rmvnorm_cholesky <- function(n, mean, cholesky){
  return(rmvnorm_cholesky_(n, mean, cholesky))
}

#'@rdname dmvnorm
#'@title Evaluate density of multivariate Normal
#'@description Function to evaluate a multivariate Normal density of each row of the given matrix.
#'@details  This function does not check anything (i.e. that the given
#' covariance is PSD). Thus it is faster than functions in standard packages, but more risky.
#'@param x is a n x d matrix of real values
#'@param mean is the mean vector
#'@param covariance is the covariance matrix
#'@return a vector of n log-density evaluations
#'@examples
#' x <- rmvnorm(10, rep(0, 5), diag(1, 5, 5))
#' dmvnorm(x, rep(-1, 5), diag(2, 5, 5))
#' # mvtnorm::dmvnorm(x, rep(-1, 5), diag(2, 5, 5), log = T)
#'@seealso \code{\link{rmvnorm}}, \code{\link{rmvnorm_cholesky}}, \code{\link{dmvnorm_cholesky_inverse}}
#'@export
dmvnorm <- function(x, mean, covariance){
  return(dmvnorm_(x, mean, covariance))
}

#'@rdname dmvnorm_cholesky_inverse
#'@title Evaluate density of multivariate Normal
#'@description Function to evaluate a multivariate Normal density of each row of the given matrix,
#' where the covariance matrix is specified through a Cholesky factor of its inverse.
#'@details  This function does not check anything (i.e. that the given
#' covariance is PSD). Thus it is faster than functions in standard packages, but more risky.
#'@param x is a n x d matrix of real values
#'@param mean is the mean vector
#'@param cholesky_inverse is a Cholesky factor of the inverse of the covariance matrix; for instance obtained
#'as \code{t(chol(solve(covariance)))}.
#'@return a vector of n log-density evaluations
#'@examples
#' x <- rmvnorm(10, rep(0, 5), diag(1, 5, 5))
#' dmvnorm_cholesky_inverse(x, rep(-1, 5), t(chol(solve(diag(2, 5, 5)))))
#' # mvtnorm::dmvnorm(x, rep(-1, 5), diag(2, 5, 5), log = T)
#'@seealso \code{\link{rmvnorm}}, \code{\link{rmvnorm_cholesky}}, \code{\link{dmvnorm}}
#'@export
dmvnorm_cholesky_inverse <- function(x, mean, cholesky_inverse){
  return(dmvnorm_cholesky_inverse_(x, mean, cholesky_inverse))
}

