source('0-random-graph-simulations.R')
set.seed(1234)
case = 'A'
n.subcases = 7
nrep = 50

for (subcase in 1:n.subcases) {
  v = seq(10, 40, by = 5)[subcase]
  K1 = K2 = 25
  K = K1 + K2
  ncomp = 2
  theta = -1.4
  alpha1 = runif(v,-0.6,0.6)
  alpha2 = runif(v,-0.6,0.6)
  for (h in 1:nrep) {
    set.seed(h)
    pop1 = vector('list', K1)
    pop2 = vector('list', K2)
    for (i in 1:K1) pop1[[i]] = random_p1(v, theta, alpha1)
    for (i in 1:K2) pop2[[i]] = random_p1(v, theta, alpha2)
    data.file = paste('data/', case, '/graphs-', case, subcase, '-', h, '.RData', sep = '')
    save.image(data.file)
  }
  print(subcase)
}
