library(igraph)

case = 'L'
n.subcases = 10
n.rep = 10

purity.avg = time.user = matrix(NA, ncol = n.subcases, nrow = n.rep)

purity.function = function(clusters, truth) {
  require(NMF)
  x = as.factor(clusters)
  out = purity(x, as.factor(truth))
  return(out)
}

make_cluster_vec = function(array, iteration) {
  K = dim(array)[1]
  prob.matrix = array[,,iteration]
  out = apply(prob.matrix, 1, function(x) which(x == max(x)))
  return(out)
}

for (aa in 1:n.subcases) {
  for (bb in 1:n.rep) {
    data.file = paste('data/', case, '/graphs-', case, aa, '-', bb, '.RData', sep = '')
    res.file = paste('results/', case, '/2-sol-', case, aa, '-', bb, '.RData', sep = '')
    load(data.file); load(res.file)
    v = seq(100, 1000, by = 100)[subcase]
    M = 2
    time.user[bb,aa] = summary(time)[1]
    logl.vec = id.max.vec = rep(NA, 10)
    clusterlist = vector('list', 10)
    for (i in 1:10) {
      #print(paste(aa, bb, i))
      id.max.vec[i] = which( sol[[i]]$loglik == max(sol[[i]]$loglik) )[1]
      logl.vec[i] = sol[[i]]$loglik[id.max.vec[i]]
    }
    x = head( which(logl.vec == max(logl.vec)), 1)
    clusters = make_cluster_vec(sol[[x]]$prob, id.max.vec[x])
    truth = c(rep(1:M, each = 25))
    purity.avg[bb,aa] = purity.function(clusters, truth)
  }
}

# median computing time
time.avg = apply(time.user, 2, median)
x = seq(100, 1000, by = 100)
lm1 = lm(time.avg ~ x+I(x^2))
plot(x, time.avg, pch = 15, xlab = 'Number of vertices (v)', ylab = 'Median computing time (seconds)')
curve(coef(lm1)[1]+coef(lm1)[2]*x+coef(lm1)[3]*x^2, add = T, lty = 3)
n.obs = x*(x-1)/2
plot(n.obs, time.avg, pch = 15, 
     xlab = 'Number of node pairs ( v(v-1)/2 )', ylab = 'Median computing time (seconds)')
lm2 = lm(time.avg ~ n.obs)
abline(a = coef(lm2)[1], b = coef(lm2)[2], lty = 3)
