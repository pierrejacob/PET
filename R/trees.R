#----------------------------------------------------------------------------------------
# Trees to store paths of particle filters
#----------------------------------------------------------------------------------------
# For more details:
# see Path storage in the particle filter, by Jacob, Murray, Rubenthaler (2015)
# [https://link.springer.com/article/10.1007/s11222-013-9445-x]
#----------------------------------------------------------------------------------------
#'@rdname tree_init
#'@title tree_init
#'@description This function initializes a tree (encoded as a \code{list}).
#'@param X_0 is a \code{dimX} by \code{N} matrix representing a set of \code{N} initial
#'\code{dimX}-dimensional particles.
#'@return This function outputs a \code{list} containing the following objects: \describe{
#'   \item{\code{N}}{Number of particles (i.e. \code{ncol(X_0)}).}
#'   \item{\code{M}}{Size of the buffer (initialized at \code{10*N}, increased automatically when needed).}
#'   \item{\code{dimX}}{Dimension of each particle (i.e. \code{nrow(X_0)}).}
#'   \item{\code{a_star}}{Vector of length \code{M} encoding ancestors (see details).}
#'   \item{\code{o_star}}{Vector of length \code{M} encoding number of offsprings (see details).}
#'   \item{\code{x_star}}{Matrix of size \code{dimX} by \code{M} encoding particles (see details).}
#'   \item{\code{l_star}}{Vector of length \code{N} containing the indices of the leaves (see details).}
#'   \item{\code{nstep}}{Height of the tree (starting at 0 with the root nodes).}
#' }
#'@details See \emph{Path storage in the particle filter}, by Jacob, Murray, Rubenthaler (2015)
#'[\url{https://link.springer.com/article/10.1007/s11222-013-9445-x}]
#'@examples
#'N = 5
#'dimX = 2
#'X_0 = matrix(rnorm(N*dimX), nrow = dimX, ncol = N)
#'tree = tree_init(X_0)
#'@export
# Initialize tree
tree_init = function(X_0){
  N = ncol(X_0)
  M = 10*N # initial buffer size M, this size is automatically increased whenever needed
  dimX = nrow(X_0)
  a_star = rep(0,M)
  o_star = rep(0,M)
  x_star = cbind(X_0,matrix(NA,nrow = dimX, ncol = M)) # dimension row-wise, particles column-wise
  l_star = (1:N)
  nstep = 0
  return (list(N = N, M = M, dimX = dimX, nstep = nstep,
               a_star = a_star, o_star = o_star, x_star = x_star, l_star = l_star))
}
#----------------------------------------------------------------------------------------
# Prune tree
# (this function is called within tree_update and does not need to be exported)
tree_prune = function(tree, o_t){
  #--- SCATTER(o_t, l_star)
  tree$o_star[tree$l_star] = o_t
  for (i in 1:tree$N){
    j = tree$l_star[i]
    while ((j>0)&&(tree$o_star[j]==0)){
      j = tree$a_star[j]
      if (j>0){
        tree$o_star[j] = tree$o_star[j] - 1
      }
    }
  }
  return (tree)
}
#----------------------------------------------------------------------------------------
# Insert new particles
# (this function is called within tree_update and does not need to be exported)
tree_insert = function(tree, X_t, a_t){
  #--- GATHER(l_star, a_t)
  b_t = tree$l_star[a_t]
  #--- TRANSFORM-PREFIX-SUM(o_star, function(x)(x==0))
  z_star = cumsum(tree$o_star==0)
  #--- LOWER-BOUND(z_star, 1:N)
  slot = 0
  i = 1
  while ((slot < tree$N) & (i < tree$M)){
    if (tree$o_star[i] == 0){
      slot = slot + 1
      tree$l_star[slot] = i
    }
    i = i+1
  } # here we can test whether slot == N, i.e. whether we have found enough slots
  if (slot < tree$N){
    # if not enough slots, double the size M
    tree$a_star = c(tree$a_star, rep(0,tree$M))
    tree$o_star = c(tree$o_star, rep(0,tree$M))
    tree$x_star = cbind(tree$x_star, matrix(0,nrow = tree$dimX, ncol = tree$M))
    for (i in (slot+1):(tree$N)){
      tree$l_star[i] = tree$M + i - slot;
    }
    tree$M = 2*(tree$M)
  }
  #--- SCATTER(b_t, l_star) and SCATTER(X_t, l_star)
  tree$a_star[tree$l_star] = b_t
  tree$x_star[,tree$l_star] = X_t
  tree$nstep = tree$nstep + 1
  return (tree)
}
#----------------------------------------------------------------------------------------
#'@rdname tree_update
#'@title tree_update
#'@description This function updates a tree given new particles and ancestors. This is achieved by
#' pruning the tree first, then inserting the new particles according to their respective ancestors.
#'@param tree is a \code{list} encoding a tree genereated via the function \code{tree_init}
#'or previously updated via the function \code{tree_update}
#'@param X_t is a \code{dimX} by \code{N} matrix representing a set of \code{N}
#'new \code{dimX}-dimensional particles.
#'@param a_t is a vector of length \code{N} containing the indices of the new particles' ancestors.
#'@return This function outputs the \code{tree} with updated fields (see details).
#'@details See documentation of \code{tree_init} for a description of the tree's fields.
#'See also \emph{Path storage in the particle filter}, by Jacob, Murray, Rubenthaler (2015)
#' [\url{https://link.springer.com/article/10.1007/s11222-013-9445-x}]
#'@examples
#'N = 5
#'dimX = 2
#'X_0 = matrix(rnorm(N*dimX), nrow = dimX, ncol = N)
#'X_1 = matrix(rnorm(N*dimX), nrow = dimX, ncol = N)
#'a_1 = sample(1:N, size = N, replace = TRUE)
#'tree = tree_init(X_0)
#'tree = tree_update(tree, X_1, a_1)
#'@export
# Update tree (prune first, then insert new particles)
tree_update = function(tree, X_t, a_t){
  tree = tree_prune(tree, o_t = tabulate(a_t, tree$N))
  tree = tree_insert(tree, X_t, a_t)
  return (tree)
}
#----------------------------------------------------------------------------------------
#'@rdname tree_getpath
#'@title tree_getpath
#'@description This function extracts a path from a tree.
#'@param tree is a \code{list} encoding a tree genereated via the function \code{tree_init} or
#'previously updated via the function \code{tree_update}.
#'@param n the index of the path to extract.
#'@return This function outputs a matrix of size \code{dimX} by \code{nstep}
#'(see documentation of \code{tree_init}) corresponding to the path indexed by \code{n}.
#'@examples
#'N = 5
#'dimX = 2
#'nstep = 10
#'
#'#--- Artificially create N particles, each of dimension dimX, evolved for nstep timesteps ---#
#'X = array(rnorm(N*dimX*(nstep+1)), dim = c(dimX,N,nstep))
#'a = replicate(nstep-1, sample(1:N, size = N, replace = TRUE))
#'tree = tree_init(X[,,1])
#'for (t in 2:nstep){
#'tree = tree_update(tree, X[,,t], a[,t-1])
#'}
#'
#'#--- List all the N final paths ---#
#'for (i in 1:N){
#'cat("Path",toString(i),"\n")
#'print(tree_getpath(tree, i))
#'}
#'@export
tree_getpath = function(tree, n){
  path = matrix(NA, nrow = tree$dimX, ncol = (tree$nstep + 1))
  j = tree$l_star[n]
  step = tree$nstep + 1
  path[,tree$nstep + 1] = tree$x_star[,j]
  step = step - 1
  while ((j >= 0)&&(step >= 0)){
    j = tree$a_star[j]
    path[,step] = tree$x_star[,j]
    step = step-1
  }
  return (path)
}

