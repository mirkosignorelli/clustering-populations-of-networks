source('0-random-graph-simulations.R')
set.seed(1234)
case = 'L'
n.subcases = 10
nrep = 10

for (subcase in 6:n.subcases) {
  nblocks = 5
  v = seq(100, 1000, by = 100)[subcase]
  K = 50
  K1 = K2 = K/2
  #
  source('0-blockprob-matrix-generator.R')
  p1 = gen.blockp.matrix(nblocks = nblocks, diag.range = c(0.2,0.5), offdiag.range = c(0.02, 0.15))
  p2 = gen.blockp.matrix(nblocks = nblocks, diag.range = c(0.2,0.5), offdiag.range = c(0.02, 0.15))
  #
  block_vector = rep(1:nblocks, each = v/nblocks)
  ncomp = 2
  for (h in 1:nrep) {
    set.seed(h)
    pop1 = vector('list', K1)
    pop2 = vector('list', K2)
    for (i in 1:K1) pop1[[i]] = random_SBM(v, block_vector, pmatrix = p1)
    for (i in 1:K2) pop2[[i]] = random_SBM(v, block_vector, pmatrix = p2)
    data.file = paste('data/', case, '/graphs-', case, subcase, '-', h, '.RData', sep = '')
    save.image(data.file)
  }
  print(subcase)
}
