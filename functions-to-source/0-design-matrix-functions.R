design_p1 = function(v, sparse = FALSE, intercept = T) {
  require('igraph')
  require('Matrix')
  a = matrix(1, v, v)
  diag(a) = 0
  b = get.edgelist(graph_from_adjacency_matrix(a, mode = 'undirected'))
  ipos = rep(1:(v*(v-1)/2),2)
  jpos = c(b[,1],b[,2])
  m = sparseMatrix(i = ipos, j = jpos, x = 1)
  trasf = sparseMatrix(i = 1:(v*(v-1)/2), j = rep(1, v*(v-1)/2), x = 1, dims = dim(m))
  trasf[,2:v] = m[,1:(v-1)] - m[,v]
  if (sparse == FALSE) trasf = as.matrix(trasf)
  if (intercept == F) trasf = trasf[,-1] # remove column with intercept
  return(trasf)
}

design_SBM = function(v, block_vector, sparse = FALSE) {
  require('igraph')
  require('Matrix')
  p = length(unique(block_vector))
  a = matrix(1, p, p)
  b = get.edgelist(graph_from_adjacency_matrix(a, mode = 'undirected'))
  out = matrix(0, v*(v-1)/2, dim(b)[1])
  k = 1
  for (i in 1:(v-1)) {
    for (j in (i+1):v) {
      position = ifelse(block_vector[i] <= block_vector[j], 
                        which(b[,1] == block_vector[i] & b[,2] == block_vector[j]),
                        which(b[,2] == block_vector[i] & b[,1] == block_vector[j]))
      out[k, position] = 1
      k = k + 1
    }
  }
  if (sparse == TRUE) out = as(out, 'sparseMatrix')
  return(out)
}

design_unconstr = function(v, sparse = FALSE) {
  require('igraph')
  require('Matrix')
  out = diag(1, v*(v-1)/2)
  if (sparse == TRUE) out = as(out, 'sparseMatrix')
  return(out)
}

# Application: functions for design matrix of degree-corrected SBM

design_dcsbm = function(v, block_vector) {
  p1.design = design_p1(v, sparse = FALSE, intercept = F)
  sbm.design = design_SBM_tdummies(v, block_vector)
  out = cbind(p1.design, sbm.design)
}

design_SBM_tdummies = function(v, block_vector, sparse = FALSE) {
  require('igraph')
  require('Matrix')
  p = length(unique(block_vector))
  out = matrix(0, v*(v-1)/2, p*(p-1)/2)
  k=1
  for (i in 1:(v-1)) {
    for (j in (i+1):v) {
      out[k,] = tail(dtrasfvec(gri = block_vector[i], grj = block_vector[j], p), p*(p-1)/2)
      k = k+1
    }
  }
  if (sparse == TRUE) out = as(out, 'sparseMatrix')
  return(out)
}

# the functions colmatr, dvec and dtrasfvec are taken from
# the scripts in:
# Signorelli, M., Wit, E. C. (2018). A penalized inference approach to stochastic 
# block modelling of community-structure in the Italian Parliament. Journal of 
# the Royal Statistical Society: Series C (Applied Statistics). 

colmatr = function(p) {
  out = matrix(NA, nrow=p, ncol = p)
  k=1
  for (i in 1:p) {
    for (j in i:p) {
      out[i,j] = k
      k = k+1
    }
  }
  return(out)
}

dvec = function(gri, grj, p) {
  dcolumn = colmatr(p)
  dvec = rep(0, p + p*(p+1)/2)
  if (gri != grj) {
    dvec[gri] = 1
    dvec[grj] = 1
  }
  else if (gri == grj) dvec[gri] = 2
  if (gri <= grj) dvec[p + dcolumn[gri,grj]] = 1
  else dvec[p + dcolumn[grj,gri]] = 1
  return(dvec)
}

dtrasfvec = function(gri, grj, p) {
  dcolumn = colmatr(p)
  dvec = dvec(gri, grj, p)
  dtrasfvec = rep(0, length = p-1 + p*(p-1)/2)
  dtrasfvec[1:(p-1)] = dvec[2:p] - dvec[1]
  k = p
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      dtrasfvec[k] = dvec[p+dcolumn[i,j]] - dvec[p+dcolumn[i,i]] - dvec[p+dcolumn[j,j]]
      k = k+1
    }
  }
  return(dtrasfvec)
}
