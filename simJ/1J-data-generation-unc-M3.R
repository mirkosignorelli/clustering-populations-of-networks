source('0-random-graph-simulations.R')
set.seed(1234)

case = 'J'
nrep = 100

v = 20
K.each = 10
ncomp = 3
K = ncomp*K.each

density = 0.3
alpha = 2
beta = alpha*(1-density)/density

for (h in 1:nrep) {
  set.seed(h)
  # define matrix of probabilities
  pvec1 = rbeta(n = v*(v-1)/2, alpha, beta)
  pvec2 = rbeta(n = v*(v-1)/2, alpha, beta)
  pvec3 = rbeta(n = v*(v-1)/2, alpha, beta)
  pmatrix1 = pmatrix2 = pmatrix3 = matrix(0, v, v)
  pmatrix1[upper.tri(pmatrix1)] = pvec1
  pmatrix2[upper.tri(pmatrix2)] = pvec2
  pmatrix3[upper.tri(pmatrix3)] = pvec3
  pmatrix1[lower.tri(pmatrix1)] = t(pmatrix1)[lower.tri(pmatrix1)]
  pmatrix2[lower.tri(pmatrix2)] = t(pmatrix2)[lower.tri(pmatrix2)]
  pmatrix3[lower.tri(pmatrix3)] = t(pmatrix3)[lower.tri(pmatrix3)]
  
  pop1 = pop2 = pop3 = vector('list', K.each)
  for (i in 1:K.each) {
    pop1[[i]] = unconstrained_model(v, pmatrix1)
    pop2[[i]] = unconstrained_model(v, pmatrix2)
    pop3[[i]] = unconstrained_model(v, pmatrix3)
  }
  data.file = paste('data/', case, '/graphs-', case, '-', h, '.RData', sep = '')
  save.image(data.file)
  cat(h)
}
