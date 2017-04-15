#'@rdname normalize_weight
#'@title Normalize vector of log-weights
#'@description Takes a vector of real values, and return the vector of normalized exponentiated values,
#' as well as average of the exponentiated but unnormalized values.
#'@export
normalize_weight <- function(logweights){
  mlw <- max(logweights)
  avew <- mlw + log(mean(exp(logweights - mlw))) # weight average
  w <- exp(logweights - mlw)
  nw <- w / sum(w) # normalized weights
  return(list(nw = nw, avew = avew))
}
