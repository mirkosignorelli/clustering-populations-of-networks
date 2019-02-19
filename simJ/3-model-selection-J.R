library(igraph)

case = 'J'
n.rep = 100
max.logl.vec = aic = bic = gic = matrix(NA, n.rep, 5)

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
  # recompute theta.hat (not saved in 2...)
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

v = 20
K.each = 10
ncomp = 3
K = ncomp*K.each
n.param = v*(v-1)/2*seq(1,5) # for unc model
n.obs = K*v*(v-1)/2

for (h in 1:n.rep) {
  # load files with optimization results
  res.file = paste('results/', case, '/2-sol-', case, '-', h, '.RData', sep = '')
  load(res.file)
  
  # get loglikelihoods
  # M = 1
  max.logl.vec[h, 1] = logLik(sol.M1)
  # M = 2
  p2 = extract.probs(sol.M2, 2)
  max.logl.vec[h, 2] = compute.loglik.undir(ylist, design, pmatrix = p2, model = 'unconstr')
  # M = 3
  p3 = extract.probs(sol.M3, 3)
  max.logl.vec[h, 3] = compute.loglik.undir(ylist, design, pmatrix = p3, model = 'unconstr')
  # M = 4
  p4 = extract.probs(sol.M4, 4)
  max.logl.vec[h, 4] = compute.loglik.undir(ylist, design, pmatrix = p4, model = 'unconstr')
  # M = 5
  p5 = extract.probs(sol.M5, 5)
  max.logl.vec[h, 5] = compute.loglik.undir(ylist, design, pmatrix = p5, model = 'unconstr')
  
  # information criteria
  aic[h, ] = -2*max.logl.vec[h, ] + 2*n.param
  bic[h, ] = -2*max.logl.vec[h, ] + log(n.obs)*n.param
  gic[h, ] = -2*max.logl.vec[h, ] + log(log(n.obs))*n.param
  if (h%%2 == 0) cat(h)
}

aic.freq = apply(aic, 1, function(x) as.numeric(which(x == min(x))))
bic.freq = apply(bic, 1, function(x) as.numeric(which(x == min(x))))
gic.freq = apply(gic, 1, function(x) as.numeric(which(x == min(x))))
table(aic.freq)
table(bic.freq)
table(gic.freq)


pmf(aic.freq, xlim = c(1,5), lwd = 3, col = 'red', title='AIC', 
    xlab = 'M', cex.main = 2, cex.axis = 1.5)
pmf(bic.freq, xlim = c(1,5), lwd = 3, col = 'red', title='BIC', 
    xlab = 'M', cex.main = 2, cex.axis = 1.5)
pmf(gic.freq, xlim = c(1,5), lwd = 3, col = 'red', title='GIC', 
    xlab = 'M', cex.main = 2, cex.axis = 1.5)

rmae = function(vec, true) mean(abs(vec-true))
rmae(aic.freq, 3)
rmae(bic.freq, 3)
rmae(gic.freq, 3)
