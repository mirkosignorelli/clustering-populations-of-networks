# NB: this script is meant to run parallelized computations
# on a cluster of computers which is based on Linux;
# the number of cores used is set here to 2, but can be changed.
# If you don't want to parallelize the computations,
# you can easily modify the code avoiding to use foreach and doMC

library(foreach)
library(doMC)
registerDoMC(cores = 2)
library(cluster)

case = 'I'
n.subcases = 6 # 2 to 7 subpopulations
n.rep = 50

for (a in 1:3) {
  for (b in 1:n.rep) {
    #
    source('0-design-matrix-functions.R')
    source('0-functions-starting-points.R')
    source('0-EM-algorithm-glm.R')
    source('0-glm-estimator.R')
    data.file = paste('data/', case, '/graphs-', case, '-', b, '.RData', sep = '')
    load(data.file)
    K = K.each*(1+a)
    ncomp = 1+a
    
    ######################################
    ###      1 : data preparation      ###
    ######################################
    adjarray = array(NA, dim = c(v, v, K))
    ylist = vector('list', length = K)
    for (i in 1:K.each) {
      adj = pop1[[i]]$adj
      adjarray[ , , i] = adj
      ylist[[i]] = adj[lower.tri(adj)]
      adj = pop2[[i]]$adj
      adjarray[ , , K.each+i] = adj
      ylist[[K.each+i]] = adj[lower.tri(adj)]
      if (a > 1) {
        adj = pop3[[i]]$adj
        adjarray[ , , 2*K.each+i] = adj
        ylist[[2*K.each+i]] = adj[lower.tri(adj)]
      }
      if (a > 2) {
        adj = pop4[[i]]$adj
        adjarray[ , , 3*K.each+i] = adj
        ylist[[3*K.each+i]] = adj[lower.tri(adj)]
      }
      if (a > 3) {
        adj = pop5[[i]]$adj
        adjarray[ , , 4*K.each+i] = adj
        ylist[[4*K.each+i]] = adj[lower.tri(adj)]
      }
      if (a > 4) {
        adj = pop6[[i]]$adj
        adjarray[ , , 5*K.each+i] = adj
        ylist[[5*K.each+i]] = adj[lower.tri(adj)]
      }
      if (a > 5) {
        adj = pop7[[i]]$adj
        adjarray[ , , 6*K.each+i] = adj
        ylist[[6*K.each+i]] = adj[lower.tri(adj)]
      }
    }
    #########################################
    ### 2: computation of starting points ###
    #########################################
    start.array = array(0, dim = c(K, ncomp, 10))
    # 1: starting point with Jaccard:
    simil = jaccard.similarity(ylist)
    km = pam(x = 1 - simil, ncomp, diss = TRUE)
    for (i in 1:ncomp) start.array[which(km$cluster==i), i, 1] = 1
    # 2: starting point with L1 distance:
    dist = l1.dist(ylist)
    km = pam(x = dist, ncomp, diss = TRUE)
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
    km = pam(x = eigendist, ncomp, diss = TRUE)
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
    if (a == 1) {
      time = system.time(
        sol <- foreach(i = 1:10) %dopar% EM_2comp(ylist, K,
                                                  start = start.array[,,i], v, block_vector = NULL,
                                                  model = 'unconstr', maxit = 100, epsilon = 1e-4)
      )      
    }
    if (a == 2) {
      time = system.time(
        sol <- foreach(i = 1:10) %dopar% EM_3comp(ylist, K,
                                                  start = start.array[,,i], v, block_vector = NULL,
                                                  model = 'unconstr', maxit = 100, epsilon = 1e-4)
      )      
    }
    if (a == 3) {
      time = system.time(
        sol <- foreach(i = 1:10) %dopar% EM_4comp(ylist, K,
                                                  start = start.array[,,i], v, block_vector = NULL,
                                                  model = 'unconstr', maxit = 100, epsilon = 1e-4)
      )      
    }
    if (a == 4) {
      time = system.time(
        sol <- foreach(i = 1:10) %dopar% EM_5comp(ylist, K,
                                                  start = start.array[,,i], v, block_vector = NULL,
                                                  model = 'unconstr', maxit = 100, epsilon = 1e-4)
      )      
    }
    if (a == 5) {
      time = system.time(
        sol <- foreach(i = 1:10) %dopar% EM_6comp(ylist, K,
                                                  start = start.array[,,i], v, block_vector = NULL,
                                                  model = 'unconstr', maxit = 100, epsilon = 1e-4)
      )      
    }
    if (a == 6) {
      time = system.time(
        sol <- foreach(i = 1:10) %dopar% EM_7comp(ylist, K,
                                                  start = start.array[,,i], v, block_vector = NULL,
                                                  model = 'unconstr', maxit = 100, epsilon = 1e-4)
      )      
    }
    save(time, sol, file = paste('results/', case, '/2-sol-', case, a, '-', b, '.RData', sep = ''))
    print(paste('Done: case', a, 'repet.', b))
  }
}
