
library(igraph)

case = 'J'
n.rep = 100
n.subcases = 4
loglik.array = array(NA, dim = c(n.rep, 5, n.subcases))

extract.probs = function(sol.list, M, niter.EM = 10) {
  EM.logl = rep(NA, niter.EM)
  for (i in 1:niter.EM) EM.logl[i] = max(sol.list[[i]]$loglik)
  best = which(EM.logl == max(EM.logl))[1]
  probs = sol.list[[best]]$prob[ , , sol.list[[best]]$niter]
  return(probs)
}

compute.loglik.undir = function(ylist, design, pmatrix, model) {
  K = dim(pmatrix)[1] 
  M = dim(pmatrix)[2]
  v = (1+sqrt(1+8*length(ylist[[1]])))/2
  n.edges = v*(v-1)/2
  # recompute theta.hat
  source('0-glm-estimator.R')
  if (M == 2) {
    mle = glm_estimator(ylist, weights = pmatrix, design, model, family = 'binomial')
    phat = cbind(mle$yhat1, mle$yhat2)
  }
  if (M == 3) {
    mle = glm_estimator_3comp(ylist, weights = pmatrix, design, model, family = 'binomial')
    phat = cbind(mle$yhat1, mle$yhat2, mle$yhat3)
  }
  if (M == 4) {
    mle = glm_estimator_4comp(ylist, weights = pmatrix, design, model, family = 'binomial')
    phat = cbind(mle$yhat1, mle$yhat2, mle$yhat3, mle$yhat4)
  }
  if (M == 5) {
    mle = glm_estimator_5comp(ylist, weights = pmatrix, design, model, family = 'binomial')
    phat = cbind(mle$yhat1, mle$yhat2, mle$yhat3, mle$yhat4, mle$yhat5)
  }
  
  y.vec = do.call(c, ylist)
  temp = rep(0, n.edges*K)
  for (m in 1:M) {
    w.vec = rep(pmatrix[ , m], each = n.edges)
    p.vec = rep(phat[ , m], K)
    ll.contr = p.vec^y.vec * (1-p.vec)^(1-y.vec)
    temp = temp + w.vec * ll.contr
  }
  out = sum(log(temp))
  return(out)
}

for (h in 1:n.rep) {
  for (l in 1:n.subcases) {
    # load files with optimization results
    res.file = paste('results/', case, '/2-sol-', case, '-', 
                     l, '-', h, '.RData', sep = '')
    temp = load(res.file)
    
    if (exists('temp')) {
      if (inherits(temp, 'try-error')) print(res.file)
      
      if (!inherits(temp, 'try-error')) {
        # get loglikelihoods
        # M = 1
        loglik.array[h, 1, l] = logLik(sol.M1)
        # M = 2
        p2 = extract.probs(sol.M2, 2)
        loglik.array[h, 2, l] = compute.loglik.undir(ylist, design, pmatrix = p2, model = 'unconstr')
        # M = 3
        p3 = extract.probs(sol.M3, 3)
        loglik.array[h, 3, l] = compute.loglik.undir(ylist, design, pmatrix = p3, model = 'unconstr')
        # M = 4
        p4 = extract.probs(sol.M4, 4)
        loglik.array[h, 4, l] = compute.loglik.undir(ylist, design, pmatrix = p4, model = 'unconstr')
        # M = 5
        p5 = extract.probs(sol.M5, 5)
        loglik.array[h, 5, l] = compute.loglik.undir(ylist, design, pmatrix = p5, model = 'unconstr')
      }
    }
        
  }
  print(h)
  if (h %% 10 == 0) save.image('results/temp-res-simJ.RData')
}

save.image('results/res-simJ.RData')



