source('0-random-graph-simulations.R')
set.seed(1234)
case = 'F'
# subcases from 1 to 6 used in script 2C, not needed here
nrep = 50

v = 30
block_vector = rep(1:3, each = v/3)
K.each = 10
# ncomp specified in script 2F

# generation of probability matrices:
pmatrix_1 = cbind(c(0.4, 0.1, 0.1), c(0.1, 0.4, 0.1), c(0.1, 0.1, 0.4))
pmatrix_2 = pmatrix_3  = pmatrix_4 = pmatrix_1
pmatrix_5 = pmatrix_6  = pmatrix_7 = pmatrix_1

pmatrix_2[1,1] = pmatrix_2[1,1] + 0.15
pmatrix_2[2,2] = pmatrix_2[2,2] - 0.15
pmatrix_3[1,1] = pmatrix_3[1,1] + 0.15
pmatrix_3[3,3] = pmatrix_3[3,3] - 0.15
pmatrix_4[2,2] = pmatrix_4[2,2] + 0.15
pmatrix_4[3,3] = pmatrix_4[3,3] - 0.15

pmatrix_5[1,1] = pmatrix_5[1,1] - 0.15
pmatrix_5[2,2] = pmatrix_5[2,2] + 0.15
pmatrix_6[1,1] = pmatrix_6[1,1] - 0.15
pmatrix_6[3,3] = pmatrix_6[3,3] + 0.15
pmatrix_7[2,2] = pmatrix_7[2,2] - 0.15
pmatrix_7[3,3] = pmatrix_7[3,3] + 0.15

for (h in 1:nrep) {
  set.seed(h)
  pop1 = pop2 = pop3 = pop4 = pop5 = pop6 = pop7 = vector('list', K.each)
  for (i in 1:K.each) {
    pop1[[i]] = random_SBM(v, block_vector, pmatrix_1)
    pop2[[i]] = random_SBM(v, block_vector, pmatrix_2)
    pop3[[i]] = random_SBM(v, block_vector, pmatrix_3)
    pop4[[i]] = random_SBM(v, block_vector, pmatrix_4)
    pop5[[i]] = random_SBM(v, block_vector, pmatrix_5)
    pop6[[i]] = random_SBM(v, block_vector, pmatrix_6)
    pop7[[i]] = random_SBM(v, block_vector, pmatrix_7)
  }
  data.file = paste('data/', case, '/graphs-', case, '-', h, '.RData', sep = '')
  save.image(data.file)
}

