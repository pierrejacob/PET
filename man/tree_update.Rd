% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trees.R
\name{tree_update}
\alias{tree_update}
\title{tree_update}
\usage{
tree_update(tree, X_t, a_t)
}
\arguments{
\item{tree}{is a \code{list} encoding a tree genereated via the function \code{tree_init}
or previously updated via the function \code{tree_update}}

\item{X_t}{is a \code{dimX} by \code{N} matrix representing a set of \code{N}
new \code{dimX}-dimensional particles.}

\item{a_t}{is a vector of length \code{N} containing the indices of the new particles' ancestors.}
}
\value{
This function outputs the \code{tree} with updated fields (see details).
}
\description{
This function updates a tree given new particles and ancestors. This is achieved by
pruning the tree first, then inserting the new particles according to their respective ancestors.
}
\details{
See documentation of \code{tree_init} for a description of the tree's fields.
See also \emph{Path storage in the particle filter}, by Jacob, Murray, Rubenthaler (2015)
[\url{https://link.springer.com/article/10.1007/s11222-013-9445-x}]
}
\examples{
N = 5
dimX = 2
X_0 = matrix(rnorm(N*dimX), nrow = dimX, ncol = N)
X_1 = matrix(rnorm(N*dimX), nrow = dimX, ncol = N)
a_1 = sample(1:N, size = N, replace = TRUE)
tree = tree_init(X_0)
tree = tree_update(tree, X_1, a_1)
}
