library(igraph)

case = 'G'
n.subcases = 7
n.rep = 50

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
    K.each = 25
    M = ncomp
    logl.vec = id.max.vec = rep(NA, 10)
    clusterlist = vector('list', 10)
    for (i in 1:10) {
      print(paste(aa, bb, i))
      id.max.vec[i] = which( sol[[i]]$loglik == max(sol[[i]]$loglik) )[1]
      logl.vec[i] = sol[[i]]$loglik[id.max.vec[i]]
    }
    x = head( which(logl.vec == max(logl.vec)), 1)
    clusters = make_cluster_vec(sol[[x]]$prob, id.max.vec[x])
    truth = c(rep(1:M, each = K.each))
    purity.avg[bb,aa] = purity.function(clusters, truth)
  }
}
    
colnames(purity.avg) = colnames(time.user) = seq(10, 40, 5)

boxplot(purity.avg, ylim = c(0.4, 1), xlab = 'v', ylab = 'purity', main = 'G', 
        cex.main = 2.3, cex.axis = 1.7, cex.lab = 1.9)
points(x = 1:7, y = rep(0.5, 7), pch = 15)
