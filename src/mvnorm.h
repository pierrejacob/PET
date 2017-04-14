#ifndef _INCL_MVNORM_
#define _INCL_MVNORM_
#include <RcppEigen.h>
using namespace Rcpp;

// generate samples from a multivariate normal distribution
NumericMatrix rmvnorm_(int nsamples, const NumericVector & mean, const NumericMatrix & covariance);
NumericMatrix rmvnorm_cholesky_(int nsamples, const NumericVector & mean, const Eigen::MatrixXd & cholesky);

// evaluate probability density function of a multivariate normal distribution
NumericVector dmvnorm_(const NumericMatrix & x, const NumericVector & mean, const NumericMatrix & covariance);
NumericVector dmvnorm_cholesky_inverse_(const NumericMatrix & x, const NumericVector & mean, const Eigen::MatrixXd & cholesky_inverse);

#endif

