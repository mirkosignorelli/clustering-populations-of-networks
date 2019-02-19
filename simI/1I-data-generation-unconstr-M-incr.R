source('0-random-graph-simulations.R')
set.seed(1234)
case = 'I'
# subcases from 1 to 6 used in script 2I, not needed here
nrep = 50

v = 15
K.each = 10
# ncomp specified in script 2I

# parameters for the beta distribution to use for \pi_{ij}:
density = 0.3
alpha = 2
beta = alpha*(1-density)/density

# define matrix of probabilities
pvec1 = rbeta(n = v*(v-1)/2, alpha, beta)
pvec2 = rbeta(n = v*(v-1)/2, alpha, beta)
pvec3 = rbeta(n = v*(v-1)/2, alpha, beta)
pvec4 = rbeta(n = v*(v-1)/2, alpha, beta)
pvec5 = rbeta(n = v*(v-1)/2, alpha, beta)
pvec6 = rbeta(n = v*(v-1)/2, alpha, beta)
pvec7 = rbeta(n = v*(v-1)/2, alpha, beta)
pmatrix1 = pmatrix2 = pmatrix3 = pmatrix4 = pmatrix5 = pmatrix6 = pmatrix7 = matrix(0, v, v)
pmatrix1[upper.tri(pmatrix1)] = pvec1
pmatrix2[upper.tri(pmatrix1)] = pvec2
pmatrix3[upper.tri(pmatrix1)] = pvec3
pmatrix4[upper.tri(pmatrix1)] = pvec4
pmatrix5[upper.tri(pmatrix1)] = pvec5
pmatrix6[upper.tri(pmatrix1)] = pvec6
pmatrix7[upper.tri(pmatrix1)] = pvec7

pmatrix1[lower.tri(pmatrix1)] = t(pmatrix1)[lower.tri(pmatrix1)]
pmatrix2[lower.tri(pmatrix2)] = t(pmatrix2)[lower.tri(pmatrix2)]
pmatrix3[lower.tri(pmatrix3)] = t(pmatrix3)[lower.tri(pmatrix3)]
pmatrix4[lower.tri(pmatrix4)] = t(pmatrix4)[lower.tri(pmatrix4)]
pmatrix5[lower.tri(pmatrix5)] = t(pmatrix5)[lower.tri(pmatrix5)]
pmatrix6[lower.tri(pmatrix6)] = t(pmatrix6)[lower.tri(pmatrix6)]
pmatrix7[lower.tri(pmatrix7)] = t(pmatrix7)[lower.tri(pmatrix7)]

for (h in 1:nrep) {
  set.seed(h)
  pop1 = pop2 = pop3 = pop4 = pop5 = pop6 = pop7 = vector('list', K.each)
  for (i in 1:K.each) {
    pop1[[i]] = unconstrained_model(v, pmatrix1)
    pop2[[i]] = unconstrained_model(v, pmatrix2)
    pop3[[i]] = unconstrained_model(v, pmatrix3)
    pop4[[i]] = unconstrained_model(v, pmatrix4)
    pop5[[i]] = unconstrained_model(v, pmatrix5)
    pop6[[i]] = unconstrained_model(v, pmatrix6)
    pop7[[i]] = unconstrained_model(v, pmatrix7)
  }
  data.file = paste('data/', case, '/graphs-', case, '-', h, '.RData', sep = '')
  save.image(data.file)
  cat(h)
}

