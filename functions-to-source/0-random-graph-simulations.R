# random_p1 generates a p1 binary graph
## alpha_constraint verifies that alpha_+ = 0 and corrects it otherwise

random_p1 = function(v, theta, alpha) {
  alpha = alpha_constraint(alpha)
  adj = matrix(0, v, v)
  for (i in 1:(v-1)) {
    for (j in (i+1):v) {
      eta = theta+alpha[i]+alpha[j]
      adj[i,j] = rbinom( 1, 1, prob = exp(eta)/(1+exp(eta)) )
    }
  }
  adj[lower.tri(adj)] = t(adj)[lower.tri(adj)]
  out = list(adj, v, theta, alpha)
  names(out) = c('adj', 'v', 'theta', 'alpha')
  return(out)
}

alpha_constraint = function(alpha) {
  total = sum(alpha)
  if (total == 0) out = alpha
  else {
    n = length(alpha)
    correction = total/n
    out = alpha - correction
  }
  return(out)
}

# random_SBM generates a SBM binary graph
random_SBM = function(v, block_vector, pmatrix) {
  adj = matrix(0, v, v)
  for (i in 1:(v-1)) {
    for (j in (i+1):v) {
      prob = pmatrix[block_vector[i], block_vector[j]]
      adj[i,j] = rbinom( 1, 1, prob = prob )
    }
  }
  adj[lower.tri(adj)] = t(adj)[lower.tri(adj)]
  out = list(adj, v, block_vector, pmatrix)
  names(out) = c('adj', 'v', 'block_vector', 'pmatrix')
  return(out)
}

# unconstrained_model generates a binary graph from a given v*v probability matrix
# design_unconstr generates the design matrix for the agnostic model 
unconstrained_model = function(v, pmatrix) {
  adj = matrix(0, v, v)
  for (i in 1:(v-1)) {
    for (j in (i+1):v) {
      adj[i,j] = rbinom( 1, 1, prob = pmatrix[i,j] )
    }
  }
  adj[lower.tri(adj)] = t(adj)[lower.tri(adj)]
  out = list(adj, v, pmatrix)
  names(out) = c('adj', 'v', 'pmatrix')
  return(out)
}

