# Jaccard index for the similarity of two vectors:
jaccard.index = function(x, y) {
  if (length(x) !=  length(y)) stop('vectors have different length')
  n11 = length(intersect(which(x == 1), which(y == 1)))
  n01 = length(intersect(which(x == 0), which(y == 1)))
  n10 = length(intersect(which(x == 1), which(y == 0)))
  out = n11 / (n11+n01+n10)
  return(out)
}

# jaccard index for a list of vectors
jaccard.similarity = function(ylist) {
  if (!is.list(ylist)) stop('ylist should be a list of vectors')
  n = length(ylist)
  if (n < 2) stop ('at least two vectors needed')
  out = matrix(NA, n, n)
  diag(out) = 1
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      out[i, j] = jaccard.index(ylist[[i]], ylist[[j]])
      out[j, i] = out[i, j]
    }
  }
  return(out)
}

# L1 distance for a list of vectors
l1.dist = function(ylist) {
  if (!is.list(ylist)) stop('ylist should be a list of vectors')
  n = length(ylist)
  if (n < 2) stop ('at least two vectors needed')
  out = matrix(NA, n, n)
  diag(out) = 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      out[i, j] = dist(rbind(ylist[[i]], ylist[[j]]), method = 'manhattan')
      out[j, i] = out[i, j]
    }
  }
  return(out)
}

# laplacian matrix (taken from FusedPCA, which is now archived on CRAN):
laplacian <- function(A, normalised = F){
  n = dim(A)[1]
  temp = apply(abs(A), 2, sum)
  D = diag(temp, nrow = n)
  temp1 = Reciprocal(sqrt(temp))
  half.D = diag(temp1, nrow = n)
  if(normalised == TRUE) 	return(half.D %*% (D - A) %*% half.D)
  if(normalised == FALSE) return(D - A)
}

Reciprocal <- function(x){
  n = length(x)
  temp = c()
  for(i in 1: n){
    if(x[i] == 0) temp[i] = 0
    else temp[i] = 1/x[i]
  }
  return(temp)
}

# reassign probabilities to some rows
reassign.prob = function(pmatrix, reass.p) {
  require('MCMCpack')
  nr = dim(pmatrix)[1]
  nc = dim(pmatrix)[2]
  n.reass = ceiling(reass.p*nr)
  rows = sample(1:nr, n.reass, replace = F)
  for (i in 1:n.reass) {
    pmatrix[rows[i], ] = rdirichlet(1, rep(1, nc))
  }
  return(pmatrix)
}

