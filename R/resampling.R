#'@export
multinomial_resampling_R <- function(nsamples, normalized_weights){
  N <- length(normalized_weights)
  return(sample(x=1:N, size=nsamples, prob=normalized_weights, replace=TRUE))
}

#'@export
multinomial_resampling <- function(nsamples, normalized_weights){
  return(multinomial_resampling_(nsamples, normalized_weights))
}

#'@export
systematic_resampling <- function(nsamples, normalized_weights){
  return(systematic_resampling_(nsamples, normalized_weights))
}
