source('0-random-graph-simulations.R')
set.seed(1234)
case = 'H'
n.subcases = 7
nrep = 50

for (subcase in 1:n.subcases) {
  v = 15
  K1 = K2 = seq(6, 30, by = 4)[subcase]
  K = K1 + K2
  ncomp = 2
  # parameters for the beta distribution to use for \pi_{ij}:
  density = 0.3
  alpha = 2
  beta = alpha*(1-density)/density
  for (h in 1:nrep) {
    set.seed(h)
    # define matrix of probabilities
    pvec1 = rbeta(n = v*(v-1)/2, alpha, beta)
    pvec2 = rbeta(n = v*(v-1)/2, alpha, beta)
    pmatrix1 = pmatrix2 = matrix(0, v, v)
    pmatrix1[upper.tri(pmatrix1)] = pvec1
    pmatrix2[upper.tri(pmatrix1)] = pvec2
    pmatrix1[lower.tri(pmatrix1)] = t(pmatrix1)[lower.tri(pmatrix1)]
    pmatrix2[lower.tri(pmatrix2)] = t(pmatrix2)[lower.tri(pmatrix2)]
    # generate the graphs
    pop1 = vector('list', K1)
    pop2 = vector('list', K2)
    for (i in 1:K1) pop1[[i]] = unconstrained_model(v, pmatrix1)
    for (i in 1:K2) pop2[[i]] = unconstrained_model(v, pmatrix2)
    data.file = paste('data/', case, '/graphs-', case, subcase, '-', h, '.RData', sep = '')
    save.image(data.file)
  }
  print(subcase)
}
