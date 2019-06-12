library(foreach)
library(doMC)
registerDoMC(cores = 10)
library(cluster)

case = 'J'
n.rep = 100
n.subcases = 4

for (h in 1:n.rep) {
  for (subcase in 1:n.subcases) {
    #
    source('0-design-matrix-functions.R')
    source('0-functions-starting-points.R')
    source('0-EM-algorithm-glm.R')
    source('0-glm-estimator.R')
    data.file = paste('data/', case, '/graphs-', case, '-',
                      subcase, '-', h, '.RData', sep = '')
    load(data.file)
    
    ######################################
    ###      1 : data preparation      ###
    ######################################
    adjarray = array(NA, dim = c(v, v, K))
    ylist = vector('list', length = K)
    for (i in 1:K.each) {
      adj = pop1[[i]]$adj
      ylist[[i]] = adj[lower.tri(adj)]
      adjarray[ , , i] = adj
      adj = pop2[[i]]$adj
      ylist[[K.each+i]] = adj[lower.tri(adj)]
      adjarray[ , , K.each+i] = adj
      adj = pop3[[i]]$adj
      ylist[[K.each*2+i]] = adj[lower.tri(adj)]
      adjarray[ , , 2*K.each+i] = adj
    }
    
    #########################################
    ### 2: computation of starting points ###
    #########################################
    start2 = array(0, dim = c(K, 2, 10))
    start3 = array(0, dim = c(K, 3, 10))
    start4 = array(0, dim = c(K, 4, 10))
    start5 = array(0, dim = c(K, 5, 10))
    # 1: starting point with Jaccard:
    simil = jaccard.similarity(ylist)
    km2 = pam(x = 1 - simil, 2, diss = TRUE)
    km3 = pam(x = 1 - simil, 3, diss = TRUE)
    km4 = pam(x = 1 - simil, 4, diss = TRUE)
    km5 = pam(x = 1 - simil, 5, diss = TRUE)
    for (i in 1:2) start2[which(km2$cluster==i), i, 1] = 1
    for (i in 1:3) start3[which(km3$cluster==i), i, 1] = 1
    for (i in 1:4) start4[which(km4$cluster==i), i, 1] = 1
    for (i in 1:5) start5[which(km5$cluster==i), i, 1] = 1
    # 2: starting point with L1 distance:
    dist = l1.dist(ylist)
    km2 = pam(x = dist, 2, diss = TRUE)
    km3 = pam(x = dist, 3, diss = TRUE)
    km4 = pam(x = dist, 4, diss = TRUE)
    km5 = pam(x = dist, 5, diss = TRUE)
    for (i in 1:2) start2[which(km2$cluster==i), i, 2] = 1
    for (i in 1:3) start3[which(km3$cluster==i), i, 2] = 1
    for (i in 1:4) start4[which(km4$cluster==i), i, 2] = 1
    for (i in 1:5) start5[which(km5$cluster==i), i, 2] = 1
    # 3: L1 distance between laplacian matrices
    eigenlapl = vector('list', length = K)
    for (i in 1:K) {
      # compute laplacian:
      lapl = laplacian(adjarray[,,i], normalised = F)
      # extract eigenvetors:
      eigenlapl[[i]] = eigen(lapl)$values
    }
    eigendist = l1.dist(eigenlapl)
    km2 = pam(x = eigendist, 2, diss = TRUE)
    km3 = pam(x = eigendist, 3, diss = TRUE)
    km4 = pam(x = eigendist, 4, diss = TRUE)
    km5 = pam(x = eigendist, 5, diss = TRUE)
    for (i in 1:2) start2[which(km2$cluster==i), i, 3] = 1
    for (i in 1:3) start3[which(km3$cluster==i), i, 3] = 1
    for (i in 1:4) start4[which(km4$cluster==i), i, 3] = 1
    for (i in 1:5) start5[which(km5$cluster==i), i, 3] = 1
    # 4-6: modification of 1
    for (i in 4:6) {
      start2[,, i] = reassign.prob(start2[,, 1], 0.3)
      start3[,, i] = reassign.prob(start3[,, 1], 0.3)
      start4[,, i] = reassign.prob(start4[,, 1], 0.3)
      start5[,, i] = reassign.prob(start5[,, 1], 0.3)
    }
    # 7-8: modification of 2
    for (i in 7:8) {
      start2[,, i] = reassign.prob(start2[,, 2], 0.3)
      start3[,, i] = reassign.prob(start3[,, 2], 0.3)
      start4[,, i] = reassign.prob(start4[,, 2], 0.3)
      start5[,, i] = reassign.prob(start5[,, 2], 0.3)
    }
    # 9-10: modification of 3
    for (i in 9:10) {
      start2[,, i] = reassign.prob(start2[,, 3], 0.3)
      start3[,, i] = reassign.prob(start3[,, 3], 0.3)
      start4[,, i] = reassign.prob(start4[,, 3], 0.3)
      start5[,, i] = reassign.prob(start5[,, 3], 0.3)
    }
    
    #######################################
    ### 3 : EM algorithm (parallelized) ###
    #######################################
    sol.M2 <- foreach(i = 1:10) %dopar% EM_2comp(ylist, K, start = start2[,,i], v, 
                                                 block_vector = NULL, model = 'unconstr', maxit = 100, epsilon = 1e-4)
    sol.M3 <- foreach(i = 1:10) %dopar% EM_3comp(ylist, K, start = start3[,,i], v, 
                                                 block_vector = NULL, model = 'unconstr', maxit = 100, epsilon = 1e-4)
    sol.M4 <- foreach(i = 1:10) %dopar% EM_4comp(ylist, K, start = start4[,,i], v, 
                                                 block_vector = NULL, model = 'unconstr', maxit = 100, epsilon = 1e-4)
    sol.M5 <- foreach(i = 1:10) %dopar% EM_5comp(ylist, K, start = start5[,,i], v, 
                                                 block_vector = NULL, model = 'unconstr', maxit = 100, epsilon = 1e-4)
    
    #compute loglik of 1 cluster only:
    yvec = do.call(c, ylist)  
    design = design_unconstr(v, sparse = FALSE)
    des = do.call( rbind, replicate(K, design, simplify=FALSE) ) 
    sol.M1 = glm(yvec ~ des-1, family = 'binomial')
    
    save(ylist, design, sol.M1, sol.M2, sol.M3, sol.M4, sol.M5, 
         file = paste('results/', case, '/2-sol-', case, '-', 
                      subcase, '-',h, '.RData', sep = ''))
    print(paste('Done: subcase', subcase, 'rep.', h))
    rm(list = setdiff(ls(), c('h', 'case', 'subcase', 'n.subcases', 'n.rep')))
    gc()
  }
}


