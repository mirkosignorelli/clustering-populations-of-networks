EM_2comp = function(ylist, K, v, start, block_vector = NULL,
                    model = c('p1', 'sbm', 'unconstr'),
                    maxit = 50, epsilon = 1e-4) {
  # block vector needed only for model == 'sbm'
  if (dim(start)[1] != K | dim(start)[2] != 2) stop('check matrix of starting values')
  if (min(start) < 0 | max(start) > 1) stop('check matrix of starting values')
  if (model == 'p1') design = design_p1(v, sparse = FALSE, intercept = F)
  else if (model == 'sbm') design = design_SBM(v, block_vector, sparse = FALSE)
  else if (model == 'unconstr') design = design_unconstr(v, sparse = FALSE)
  ncomp = 2
  membership_prob = array(dim = c(K, ncomp, maxit))
  loglik = rep(NA, maxit)
  i = 1
  logl_change = Inf
  while (logl_change > epsilon & i <= maxit) {
    # E step
    if (i == 1) {
      membership_prob[,,i] = start
    }
    else if (i > 1) {
      graphprob1 = numeric();graphprob2 = numeric()
      for (j in 1:K) {
        obs = ylist[[j]]
        e = v*(v-1)/2
        lprob1 = dbinom(obs, 1, glm.est$yhat1, log = T)
        graphprob1[j] = sum(lprob1)
        lprob2 = dbinom(obs, 1, glm.est$yhat2, log = T)
        graphprob2[j] = sum(lprob2)
      }
      log.lik.matrix = cbind(graphprob1, graphprob2)
      membership_prob[,,i] = prob_Kcomp(log.lik.matrix)
    }
    # M step
    current.weights = membership_prob[,,i]
    glm.est = glm_estimator(ylist = ylist, weights = current.weights, design = design, 
                            model = model, family = 'binomial')
    loglik[i] = glm.est$loglik
    if (i > 1) logl_change = loglik[i] - loglik[i-1]
    # next iteration
    i = i + 1
  }
  niter = i-1
  out = list('niter' = niter, 'prob' = membership_prob[,,1:niter], 'loglik' = loglik[1:niter])
  return(out)
}

EM_3comp = function(ylist, K, v, start, block_vector = NULL,
                    model = c('p1', 'sbm', 'unconstr'),
                    maxit = 50, epsilon = 1e-4) {
  # block vector needed only for model == 'sbm'
  if (dim(start)[1] != K | dim(start)[2] != 3) stop('check matrix of starting values')
  if (min(start) < 0 | max(start) > 1) stop('check matrix of starting values')
  if (model == 'p1') design = design_p1(v, sparse = FALSE, intercept = F)
  else if (model == 'sbm') design = design_SBM(v, block_vector, sparse = FALSE)
  else if (model == 'unconstr') design = design_unconstr(v, sparse = FALSE)
  ncomp = 3
  membership_prob = array(dim = c(K, ncomp, maxit))
  loglik = rep(NA, maxit)
  i = 1
  logl_change = Inf
  while (logl_change > epsilon & i <= maxit) {
    # E step
    if (i == 1) {
      membership_prob[,,i] = start
    }
    else if (i > 1) {
      graphprob1 = graphprob2 = graphprob3 = numeric()
      for (j in 1:K) {
        obs = ylist[[j]]
        e = v*(v-1)/2
        lprob1 = dbinom(obs, 1, glm.est$yhat1, log = T)
        graphprob1[j] = sum(lprob1)
        lprob2 = dbinom(obs, 1, glm.est$yhat2, log = T)
        graphprob2[j] = sum(lprob2)
        lprob3 = dbinom(obs, 1, glm.est$yhat3, log = T)
        graphprob3[j] = sum(lprob3)
      }
      log.lik.matrix = cbind(graphprob1, graphprob2, graphprob3)
      membership_prob[,,i] = prob_Kcomp(log.lik.matrix)
    }
    # M step
    current.weights = membership_prob[,,i]
    glm.est = glm_estimator_3comp(ylist = ylist, weights = current.weights, design = design, 
                                  model = model, family = 'binomial')
    loglik[i] = glm.est$loglik
    if (i > 1) logl_change = loglik[i] - loglik[i-1]
    # next iteration
    i = i + 1
  }
  niter = i-1
  out = list('niter' = niter, 'prob' = membership_prob[,,1:niter], 'loglik' = loglik[1:niter])
  return(out)
}

EM_4comp = function(ylist, K, v, start, block_vector = NULL,
                    model = c('p1', 'sbm', 'unconstr'),
                    maxit = 50, epsilon = 1e-4) {
  if (dim(start)[1] != K | dim(start)[2] != 4) stop('check matrix of starting values')
  if (min(start) < 0 | max(start) > 1) stop('check matrix of starting values')
  # block vector needed only for model == 'sbm'
  if (model == 'p1') design = design_p1(v, sparse = FALSE, intercept = F)
  else if (model == 'sbm') design = design_SBM(v, block_vector, sparse = FALSE)
  else if (model == 'unconstr') design = design_unconstr(v, sparse = FALSE)
  ncomp = 4
  membership_prob = array(dim = c(K, ncomp, maxit))
  loglik = rep(NA, maxit)
  i = 1
  logl_change = Inf
  while (logl_change > epsilon & i <= maxit) {
    # E step
    if (i == 1) {
      membership_prob[,,i] = start
    }
    else if (i > 1) {
      graphprob1 = graphprob2 = graphprob3 = graphprob4 = numeric()
      for (j in 1:K) {
        obs = ylist[[j]]
        e = v*(v-1)/2
        lprob1 = dbinom(obs, 1, glm.est$yhat1, log = T)
        graphprob1[j] = sum(lprob1)
        lprob2 = dbinom(obs, 1, glm.est$yhat2, log = T)
        graphprob2[j] = sum(lprob2)
        lprob3 = dbinom(obs, 1, glm.est$yhat3, log = T)
        graphprob3[j] = sum(lprob3)
        lprob4 = dbinom(obs, 1, glm.est$yhat4, log = T)
        graphprob4[j] = sum(lprob4)
      }
      log.lik.matrix = cbind(graphprob1, graphprob2, graphprob3, graphprob4)
      membership_prob[,,i] = prob_Kcomp(log.lik.matrix)
    }
    # M step
    current.weights = membership_prob[,,i]
    glm.est = glm_estimator_4comp(ylist = ylist, weights = current.weights, design = design, 
                                  model = model, family = 'binomial')
    loglik[i] = glm.est$loglik
    if (i > 1) logl_change = loglik[i] - loglik[i-1]
    # next iteration
    i = i + 1
  }
  niter = i-1
  out = list('niter' = niter, 'prob' = membership_prob[,,1:niter], 'loglik' = loglik[1:niter])
  return(out)
}

EM_5comp = function(ylist, K, v, start, block_vector = NULL,
                    model = c('p1', 'sbm', 'unconstr'),
                    maxit = 50, epsilon = 1e-4) {
  # block vector needed only for model == 'sbm'
  if (dim(start)[1] != K | dim(start)[2] != 5) stop('check matrix of starting values')
  if (min(start) < 0 | max(start) > 1) stop('check matrix of starting values')
  if (model == 'p1') design = design_p1(v, sparse = FALSE, intercept = F)
  else if (model == 'sbm') design = design_SBM(v, block_vector, sparse = FALSE)
  else if (model == 'unconstr') design = design_unconstr(v, sparse = FALSE)
  ncomp = 5
  membership_prob = array(dim = c(K, ncomp, maxit))
  loglik = rep(NA, maxit)
  i = 1
  logl_change = Inf
  while (logl_change > epsilon & i <= maxit) {
    # E step
    if (i == 1) {
      membership_prob[,,i] = start
    }
    else if (i > 1) {
      graphprob1 = graphprob2 = graphprob3 = graphprob4 = graphprob5 = numeric()
      for (j in 1:K) {
        obs = ylist[[j]]
        e = v*(v-1)/2
        lprob1 = dbinom(obs, 1, glm.est$yhat1, log = T)
        graphprob1[j] = sum(lprob1)
        lprob2 = dbinom(obs, 1, glm.est$yhat2, log = T)
        graphprob2[j] = sum(lprob2)
        lprob3 = dbinom(obs, 1, glm.est$yhat3, log = T)
        graphprob3[j] = sum(lprob3)
        lprob4 = dbinom(obs, 1, glm.est$yhat4, log = T)
        graphprob4[j] = sum(lprob4)
        lprob5 = dbinom(obs, 1, glm.est$yhat5, log = T)
        graphprob5[j] = sum(lprob5)
      }
      log.lik.matrix = cbind(graphprob1, graphprob2, graphprob3, graphprob4, graphprob5)
      membership_prob[,,i] = prob_Kcomp(log.lik.matrix)
    }
    # M step
    current.weights = membership_prob[,,i]
    glm.est = glm_estimator_5comp(ylist = ylist, weights = current.weights, design = design, 
                                  model = model, family = 'binomial')
    loglik[i] = glm.est$loglik
    if (i > 1) logl_change = loglik[i] - loglik[i-1]
    # next iteration
    i = i + 1
  }
  niter = i-1
  out = list('niter' = niter, 'prob' = membership_prob[,,1:niter], 'loglik' = loglik[1:niter])
  return(out)
}

EM_6comp = function(ylist, K, v, start, block_vector = NULL,
                    model = c('p1', 'sbm', 'unconstr'),
                    maxit = 50, epsilon = 1e-4) {
  # block vector needed only for model == 'sbm'
  if (dim(start)[1] != K | dim(start)[2] != 6) stop('check matrix of starting values')
  if (min(start) < 0 | max(start) > 1) stop('check matrix of starting values')
  if (model == 'p1') design = design_p1(v, sparse = FALSE, intercept = F)
  else if (model == 'sbm') design = design_SBM(v, block_vector, sparse = FALSE)
  else if (model == 'unconstr') design = design_unconstr(v, sparse = FALSE)
  ncomp = 6
  membership_prob = array(dim = c(K, ncomp, maxit))
  loglik = rep(NA, maxit)
  i = 1
  logl_change = Inf
  while (logl_change > epsilon & i <= maxit) {
    # E step
    if (i == 1) {
      membership_prob[,,i] = start
    }
    else if (i > 1) {
      graphprob1 = graphprob2 = graphprob3 = graphprob4 = graphprob5 = graphprob6 = numeric()
      for (j in 1:K) {
        obs = ylist[[j]]
        e = v*(v-1)/2
        lprob1 = dbinom(obs, 1, glm.est$yhat1, log = T)
        graphprob1[j] = sum(lprob1)
        lprob2 = dbinom(obs, 1, glm.est$yhat2, log = T)
        graphprob2[j] = sum(lprob2)
        lprob3 = dbinom(obs, 1, glm.est$yhat3, log = T)
        graphprob3[j] = sum(lprob3)
        lprob4 = dbinom(obs, 1, glm.est$yhat4, log = T)
        graphprob4[j] = sum(lprob4)
        lprob5 = dbinom(obs, 1, glm.est$yhat5, log = T)
        graphprob5[j] = sum(lprob5)
        lprob6 = dbinom(obs, 1, glm.est$yhat6, log = T)
        graphprob6[j] = sum(lprob6)
      }
      log.lik.matrix = cbind(graphprob1, graphprob2, graphprob3, graphprob4, graphprob5, graphprob6)
      membership_prob[,,i] = prob_Kcomp(log.lik.matrix)
    }
    # M step
    current.weights = membership_prob[,,i]
    glm.est = glm_estimator_6comp(ylist = ylist, weights = current.weights, design = design, 
                                  model = model, family = 'binomial')
    loglik[i] = glm.est$loglik
    if (i > 1) logl_change = loglik[i] - loglik[i-1]
    # next iteration
    i = i + 1
  }
  niter = i-1
  out = list('niter' = niter, 'prob' = membership_prob[,,1:niter], 'loglik' = loglik[1:niter])
  return(out)
}

EM_7comp = function(ylist, K, v, start, block_vector = NULL,
                    model = c('p1', 'sbm', 'unconstr'),
                    maxit = 50, epsilon = 1e-4) {
  # block vector needed only for model == 'sbm'
  if (dim(start)[1] != K | dim(start)[2] != 7) stop('check matrix of starting values')
  if (min(start) < 0 | max(start) > 1) stop('check matrix of starting values')
  if (model == 'p1') design = design_p1(v, sparse = FALSE, intercept = F)
  else if (model == 'sbm') design = design_SBM(v, block_vector, sparse = FALSE)
  else if (model == 'unconstr') design = design_unconstr(v, sparse = FALSE)
  ncomp = 7
  membership_prob = array(dim = c(K, ncomp, maxit))
  loglik = rep(NA, maxit)
  i = 1
  logl_change = Inf
  while (logl_change > epsilon & i <= maxit) {
    # E step
    if (i == 1) {
      membership_prob[,,i] = start
    }
    else if (i > 1) {
      graphprob1 = graphprob2 = graphprob3 = graphprob4 = graphprob5 = graphprob6 = graphprob7 = numeric()
      for (j in 1:K) {
        obs = ylist[[j]]
        e = v*(v-1)/2
        lprob1 = dbinom(obs, 1, glm.est$yhat1, log = T)
        graphprob1[j] = sum(lprob1)
        lprob2 = dbinom(obs, 1, glm.est$yhat2, log = T)
        graphprob2[j] = sum(lprob2)
        lprob3 = dbinom(obs, 1, glm.est$yhat3, log = T)
        graphprob3[j] = sum(lprob3)
        lprob4 = dbinom(obs, 1, glm.est$yhat4, log = T)
        graphprob4[j] = sum(lprob4)
        lprob5 = dbinom(obs, 1, glm.est$yhat5, log = T)
        graphprob5[j] = sum(lprob5)
        lprob6 = dbinom(obs, 1, glm.est$yhat6, log = T)
        graphprob6[j] = sum(lprob6)
        lprob7 = dbinom(obs, 1, glm.est$yhat7, log = T)
        graphprob7[j] = sum(lprob7)
      }
      log.lik.matrix = cbind(graphprob1, graphprob2, graphprob3, graphprob4, graphprob5, graphprob6, graphprob7)
      membership_prob[,,i] = prob_Kcomp(log.lik.matrix)
    }
    # M step
    current.weights = membership_prob[,,i]
    glm.est = glm_estimator_7comp(ylist = ylist, weights = current.weights, design = design, 
                                  model = model, family = 'binomial')
    loglik[i] = glm.est$loglik
    if (i > 1) logl_change = loglik[i] - loglik[i-1]
    # next iteration
    i = i + 1
  }
  niter = i-1
  out = list('niter' = niter, 'prob' = membership_prob[,,1:niter], 'loglik' = loglik[1:niter])
  return(out)
}

prob_Kcomp = function(log.p.graphs) {
  q = dim(log.p.graphs)
  K = q[1]
  M = q[2]
  out = matrix(0, nrow = K, ncol = M)
  for (i in 1:K) {
    nums = exp(log.p.graphs[i,])
    den = sum(nums)
    if (den == Inf | is.nan(den)) {
      which.inf = which(is.infinite(nums) | is.nan(nums))
      if (length(which.inf) > 1) stop('two clusters with infinite numerators for same graph')
      else out[i,which.inf] = 1
    }
    else if (den == 0) {
      temp = which(log.p.graphs[i,] == max(log.p.graphs[i,]))
      out[i, temp] = 1
    }
    else out[i,] = nums / den
  }
  return(out)  
}
