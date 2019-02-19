source('0-random-graph-simulations.R')
set.seed(1234)
case = 'D'
n.subcases = 7
nrep = 50

pmatrix_1 = cbind(c(0.4, 0.1, 0.1), c(0.1, 0.5, 0.2), c(0.1, 0.2, 0.3))
pmatrix_2 = cbind(c(0.3, 0.2, 0.2), c(0.2, 0.3, 0.1), c(0.2, 0.1, 0.4))

for (subcase in 1:n.subcases) {
  v = seq(10, 40, by = 5)[subcase]
  block_vector = rep(1:3, each = v/3)
  if (v - length(block_vector) ==1) block_vector = c(block_vector, 1)
  if (v - length(block_vector) ==2) block_vector = c(block_vector, 1, 2)
  K1 = K2 = 25
  K = K1 + K2
  ncomp = 2
  for (h in 1:nrep) {
    set.seed(h)
    pop1 = vector('list', K1)
    pop2 = vector('list', K2)
    for (i in 1:K1) pop1[[i]] = random_SBM(v, block_vector, pmatrix = pmatrix_1)
    for (i in 1:K2) pop2[[i]] = random_SBM(v, block_vector, pmatrix = pmatrix_2)
    data.file = paste('data/', case, '/graphs-', case, subcase, '-', h, '.RData', sep = '')
    save.image(data.file)
  }
  print(subcase)
}
