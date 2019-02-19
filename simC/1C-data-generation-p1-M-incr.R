source('0-random-graph-simulations.R')
set.seed(1234)
case = 'C'
# subcases from 1 to 6 used in script 2C, not needed here
nrep = 50

v = 30
K.each = 10
# ncomp specified in script 2C

theta = -1.4
alpha1 = runif(v,-0.6,0.6)
alpha2 = runif(v,-0.6,0.6)
alpha3 = runif(v,-0.6,0.6)
alpha4 = runif(v,-0.6,0.6)
alpha5 = runif(v,-0.6,0.6)
alpha6 = runif(v,-0.6,0.6)
alpha7 = runif(v,-0.6,0.6)

for (h in 1:nrep) {
  set.seed(h)
  pop1 = pop2 = pop3 = pop4 = pop5 = pop6 = pop7 = vector('list', K.each)
  for (i in 1:K.each) {
    pop1[[i]] = random_p1(v, theta, alpha1)
    pop2[[i]] = random_p1(v, theta, alpha2)
    pop3[[i]] = random_p1(v, theta, alpha3)
    pop4[[i]] = random_p1(v, theta, alpha4)
    pop5[[i]] = random_p1(v, theta, alpha5)
    pop6[[i]] = random_p1(v, theta, alpha6)
    pop7[[i]] = random_p1(v, theta, alpha7)
  }
  data.file = paste('data/', case, '/graphs-', case, '-', h, '.RData', sep = '')
  save.image(data.file)
  cat(h)
}

