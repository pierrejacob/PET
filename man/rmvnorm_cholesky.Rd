% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivariatenormal.R
\name{rmvnorm_cholesky}
\alias{rmvnorm_cholesky}
\title{Generate multivariate Normals}
\usage{
rmvnorm_cholesky(n, mean, cholesky)
}
\arguments{
\item{n}{is the number of desired samples}

\item{mean}{is the mean vector}

\item{cholesky}{is the cholesky factor of the covariance matrix, obtained e.g. with \code{chol()}.}
}
\value{
a n x d matrix where d is both the length of \code{mean} and the dimension of \code{cholesky}. Each row contains a draw from
the desired multivariate Normal.
}
\description{
Function to generate draws from a multivariate Normal distribution, given Cholesky factor of the covariance.
}
\details{
This function does not check anything (i.e. that the given
covariance is PSD). Thus it is faster than functions in standard packages, but more risky.
}
\examples{
rmvnorm_cholesky(10, rep(0, 5), chol(diag(1, 5, 5)))
}
\seealso{
\code{\link{rmvnorm}}, \code{\link{dmvnorm}}, \code{\link{dmvnorm_cholesky_inverse}}
}
