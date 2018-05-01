#'@rdname multinomial_resampling_R
#'@title Multinomial resampling (slow implementation)
#'@description Sample \code{nsamples} times from a categorical
#' distribution with probabilities given by \code{normalized_weights}.
#'@param nsamples is a number of desired samples
#'@param normalized_weights is a vector of probabilities (non-negative values summing to one)
#'@return a vector of \code{nsamples} integers between 1 and the size of \code{normalized_weights}.
#'@examples
#' N <- 1000
#' logweights <- rnorm(N)
#' normalize_weight_results <- normalize_weight(logweights)
#' normalized_weights <- normalize_weight_results$nw
#' multinomial_resampling_R(10, normalized_weights)
#'@export
multinomial_resampling_R <- function(nsamples, normalized_weights){
  N <- length(normalized_weights)
  return(sample(x=1:N, size=nsamples, prob=normalized_weights, replace=TRUE))
}

#'@rdname multinomial_resampling
#'@title Multinomial resampling (fast implementation)
#'@description Sample \code{nsamples} times from a categorical
#' distribution with probabilities given by \code{normalized_weights}.
#'@details This implementation is in C++, based on drawing sorted uniforms.
#'@param nsamples is a number of desired samples
#'@param normalized_weights is a vector of probabilities (non-negative values summing to one)
#'@return a vector of \code{nsamples} integers between 1 and the size of \code{normalized_weights}.
#'@examples
#' N <- 1000
#' logweights <- rnorm(N)
#' normalize_weight_results <- normalize_weight(logweights)
#' normalized_weights <- normalize_weight_results$nw
#' multinomial_resampling(10, normalized_weights)
#'@export
multinomial_resampling <- function(nsamples, normalized_weights){
  return(multinomial_resampling_(nsamples, normalized_weights))
}


#'@rdname systematic_resampling
#'@title Systematic resampling
#'@description Systematic resampling returns a vector of ancestors based on the given
#'weights, using only one uniform.
#'@details See Comparison of Resampling Schemes for Particle Filtering, Douc, Cappe and Moulines,
#' https://arxiv.org/pdf/cs/0507025.pdf
#'@param nsamples is a number of desired samples
#'@param normalized_weights is a vector of probabilities (non-negative values summing to one)
#'@return a vector of \code{nsamples} integers between 1 and the size of \code{normalized_weights}.
#'@examples
#' N <- 1000
#' logweights <- rnorm(N)
#' normalize_weight_results <- normalize_weight(logweights)
#' normalized_weights <- normalize_weight_results$nw
#' systematic_resampling(10, normalized_weights)
#'@export
systematic_resampling <- function(nsamples, normalized_weights){
  return(systematic_resampling_(nsamples, normalized_weights))
}


#'@rdname ssp_resampling
#'@title SSP resampling
#'@description  Srinivasan Sampling Process (SSP) resampling returns a vector of \code{nsamples} ancestors
#'(integers from 1 to \code{N}) based on a given vector of \code{N} normalized weights.
#'This function internally uses a vector of \code{N} independent Unif(0,1).
#'
#'The SSP resampling is unbiased: if \code{A = ssp_resampling(N,W)}, then
#'for every \code{i = 1,...,N}, the expectation of \code{sum(A==i)} is equal to \code{N*W[i]}.
#'
#'The SSP resampling has the additional property that \code{sum(A==i)} is exactly equal to either
#' \code{floor(N*W[i])} or \code{(1+floor(N*W[i]))}.
#'@details See Negative association, ordering and convergence of resampling methods, by Gerber,
#'Chopin, and Whiteley (2017) [https://arxiv.org/abs/1707.01845]
#'@param normalized_weights is a vector of probabilities (non-negative values summing to one)
#'@return A vector of ancestors of the same length as \code{normalized_weights}
#'@examples
#' N <- 1000
#' logweights <- rnorm(N)
#' normalize_weight_results <- normalize_weight(logweights)
#' normalized_weights <- normalize_weight_results$nw
#' ssp_resampling(normalized_weights)
#'@export
ssp_resampling <- function(nsamples, normalized_weights){
  return(SSP_resampling_(nsamples, normalized_weights))
}
