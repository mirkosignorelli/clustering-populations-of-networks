# NB: this script is meant to run parallelized computations
# on a cluster of computers which is based on Linux;
# the number of cores used is set here to 2, but can be changed.
# If you don't want to parallelize the computations,
# you can easily modify the code avoiding to use foreach and doMC

library(foreach)
library(doMC)
registerDoMC(cores = 2)
library(cluster)

case = 'K'
n.subcases = 10
n.rep = 10

for (a in 1:n.subcases) {
  for (b in 1:n.rep) {
#
    source('0-design-matrix-functions.R')
    source('0-functions-starting-points.R')
    source('0-EM-algorithm-glm.R')
    source('0-glm-estimator.R')
    data.file = paste('data/', case, '/graphs-', case, a, '-', b, '.RData', sep = '')
    load(data.file)
    
    ######################################
    ###      1 : data preparation      ###
    ######################################
    adjarray = array(NA, dim = c(v, v, K))
    ylist = vector('list', length = K)
    for (i in 1:K1) {
      adj = pop1[[i]]$adj
      adjarray[ , , i] = adj
      ylist[[i]] = adj[lower.tri(adj)]
    }
    for (i in 1:K2) {
      adj = pop2[[i]]$adj
      adjarray[ , , K1+i] = adj
      ylist[[K1+i]] = adj[lower.tri(adj)]
    }
    
    #########################################
    ### 2: computation of starting points ###
    #########################################
    start.array = array(0, dim = c(K, ncomp, 10))
    # 1: starting point with Jaccard:
    simil = jaccard.similarity(ylist)
    km = pam(x = 1 - simil, 2, diss = TRUE)
    for (i in 1:ncomp) start.array[which(km$cluster==i), i, 1] = 1
    # 2: starting point with L1 distance:
    dist = l1.dist(ylist)
    km = pam(x = dist, 2, diss = TRUE)
    for (i in 1:ncomp) start.array[which(km$cluster==i), i, 2] = 1
    # 3: L1 distance between laplacian matrices
    eigenlapl = vector('list', length = K)
    for (i in 1:K) {
      # compute laplacian:
      lapl = laplacian(adjarray[,,i], normalised = F)
      # extract eigenvetors:
      eigenlapl[[i]] = eigen(lapl)$values
    }
    eigendist = l1.dist(eigenlapl)
    km = pam(x = eigendist, 2, diss = TRUE)
    for (i in 1:ncomp) start.array[which(km$cluster==i), i, 3] = 1
    # 4-6: modification of 1
    for (i in 4:6) start.array[,, i] = reassign.prob(start.array[,, 1], 0.3)
    # 7-8: modification of 2
    for (i in 7:8) start.array[,, i] = reassign.prob(start.array[,, 2], 0.3)
    # 9-10: modification of 3
    for (i in 9:10) start.array[,, i] = reassign.prob(start.array[,, 3], 0.3)
    
    #######################################
    ### 3 : EM algorithm (parallelized) ###
    #######################################
    # EM algorithm (parallelized):
    time = system.time(
      sol <- foreach(i = 1:10) %dopar% EM_2comp(ylist, K,
                  start = start.array[,,i], v, block_vector = block_vector,
                  model = 'sbm', maxit = 50, epsilon = 1e-4)
    )
    
    save(time, sol, file = paste('results/', case, '/2-sol-', case, a, '-', b, '.RData', sep = ''))
    print(paste('Done: case', a, 'repet.', b))
    rm(list = setdiff(ls(), c('a', 'b', 'case', 'n.subcases', 'n.rep')))    
  }
}



