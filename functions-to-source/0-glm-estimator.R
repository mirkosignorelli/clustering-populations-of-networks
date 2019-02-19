# glm estimator with 2 components, p1 or sbm models:
glm_estimator = function(ylist, weights, design, model = c('p1', 'sbm', 'unconstr'), family = 'binomial') {
  if (!exists('v')) v = (1+sqrt(1+8*length(ylist[[1]])))/2
  w1 = rep(weights[,1], each = v*(v-1)/2)
  w2 = rep(weights[,2], each = v*(v-1)/2)
  y = do.call(c, ylist)
  K = length(ylist) # number of graphs
  X = do.call( rbind, replicate(K, design, simplify=FALSE) )  

  if (model == 'p1') {
    comp1 = glm(y ~ X, weights = w1, family = family)
    comp2 = glm(y ~ X, weights = w2, family = family)
  }
  else if (model %in% c('sbm','unconstr')) {
    comp1 = glm(y ~ X-1, weights = w1, family = family)
    comp2 = glm(y ~ X-1, weights = w2, family = family)
  } 

  if (model == 'p1') newdata = cbind(1, design)
  else if (model %in% c('sbm','unconstr')) newdata = design
  yhat1 = plogis(newdata%*%coef(comp1))
  yhat2 = plogis(newdata%*%coef(comp2))
  loglik = logLik(comp1) + logLik(comp2)
  return(list('comp1'=comp1, 'comp2'=comp2, 'yhat1' = yhat1, 'yhat2' = yhat2, 'loglik' = as.numeric(loglik)))
}

# glm estimator with 3 components, p1 or sbm models:
glm_estimator_3comp = function(ylist, weights, design, model = c('p1', 'sbm', 'unconstr'), family = 'binomial') {
  if (!exists('v')) v = (1+sqrt(1+8*length(ylist[[1]])))/2
  w1 = rep(weights[,1], each = v*(v-1)/2)
  w2 = rep(weights[,2], each = v*(v-1)/2)
  w3 = rep(weights[,3], each = v*(v-1)/2)
  y = do.call(c, ylist)
  K = length(ylist) # number of graphs
  X = do.call( rbind, replicate(K, design, simplify=FALSE) )  
  
  if (model == 'p1') {
    comp1 = glm(y ~ X, weights = w1, family = family)
    comp2 = glm(y ~ X, weights = w2, family = family)
    comp3 = glm(y ~ X, weights = w3, family = family)
  }
  else if (model %in% c('sbm','unconstr')) {
    comp1 = glm(y ~ X-1, weights = w1, family = family)
    comp2 = glm(y ~ X-1, weights = w2, family = family)
    comp3 = glm(y ~ X-1, weights = w3, family = family)
  }   
    
  if (model == 'p1') newdata = cbind(1, design)
  else if (model %in% c('sbm','unconstr')) newdata = design
  yhat1 = plogis(newdata%*%coef(comp1))
  yhat2 = plogis(newdata%*%coef(comp2))
  yhat3 = plogis(newdata%*%coef(comp3))
  loglik = logLik(comp1)+logLik(comp2)+logLik(comp3)
  return(list('comp1'=comp1, 'comp2'=comp2, 'comp3'=comp3,
              'yhat1' = yhat1, 'yhat2' = yhat2, 'yhat3' = yhat3,
              'loglik' = as.numeric(loglik)))
}

# glm estimator with 4 components, p1 or sbm models:
glm_estimator_4comp = function(ylist, weights, design, model = c('p1', 'sbm', 'unconstr'), family = 'binomial') {
  if (!exists('v')) v = (1+sqrt(1+8*length(ylist[[1]])))/2
  w1 = rep(weights[,1], each = v*(v-1)/2)
  w2 = rep(weights[,2], each = v*(v-1)/2)
  w3 = rep(weights[,3], each = v*(v-1)/2)
  w4 = rep(weights[,4], each = v*(v-1)/2)
  y = do.call(c, ylist)
  K = length(ylist) # number of graphs
  X = do.call( rbind, replicate(K, design, simplify=FALSE) )  
  
  if (model == 'p1') {
    comp1 = glm(y ~ X, weights = w1, family = family)
    comp2 = glm(y ~ X, weights = w2, family = family)
    comp3 = glm(y ~ X, weights = w3, family = family)
    comp4 = glm(y ~ X, weights = w4, family = family)
  }
  else if (model %in% c('sbm','unconstr')) {
    comp1 = glm(y ~ X-1, weights = w1, family = family)
    comp2 = glm(y ~ X-1, weights = w2, family = family)
    comp3 = glm(y ~ X-1, weights = w3, family = family)
    comp4 = glm(y ~ X-1, weights = w4, family = family)
  }   
  
  if (model == 'p1') newdata = cbind(1, design)
  else if (model %in% c('sbm','unconstr')) newdata = design
  yhat1 = plogis(newdata%*%coef(comp1))
  yhat2 = plogis(newdata%*%coef(comp2))
  yhat3 = plogis(newdata%*%coef(comp3))
  yhat4 = plogis(newdata%*%coef(comp4))
  loglik = logLik(comp1)+logLik(comp2)+logLik(comp3)+logLik(comp4)
  return(list('comp1'=comp1, 'comp2'=comp2, 'comp3'=comp3, 'comp4'=comp4,
              'yhat1' = yhat1, 'yhat2' = yhat2, 'yhat3' = yhat3, 'yhat4' = yhat4, 
              'loglik' = as.numeric(loglik)))
}

# glm estimator with 5 components, p1 or sbm models:
glm_estimator_5comp = function(ylist, weights, design, model = c('p1', 'sbm', 'unconstr'), family = 'binomial') {
  if (!exists('v')) v = (1+sqrt(1+8*length(ylist[[1]])))/2
  w1 = rep(weights[,1], each = v*(v-1)/2)
  w2 = rep(weights[,2], each = v*(v-1)/2)
  w3 = rep(weights[,3], each = v*(v-1)/2)
  w4 = rep(weights[,4], each = v*(v-1)/2)
  w5 = rep(weights[,5], each = v*(v-1)/2)
  y = do.call(c, ylist)
  K = length(ylist) # number of graphs
  X = do.call( rbind, replicate(K, design, simplify=FALSE) )  
  
  if (model == 'p1') {
    comp1 = glm(y ~ X, weights = w1, family = family)
    comp2 = glm(y ~ X, weights = w2, family = family)
    comp3 = glm(y ~ X, weights = w3, family = family)
    comp4 = glm(y ~ X, weights = w4, family = family)
    comp5 = glm(y ~ X, weights = w5, family = family)
  }
  else if (model %in% c('sbm','unconstr')) {
    comp1 = glm(y ~ X-1, weights = w1, family = family)
    comp2 = glm(y ~ X-1, weights = w2, family = family)
    comp3 = glm(y ~ X-1, weights = w3, family = family)
    comp4 = glm(y ~ X-1, weights = w4, family = family)
    comp5 = glm(y ~ X-1, weights = w5, family = family)
  }   
  
  if (model == 'p1') newdata = cbind(1, design)
  else if (model %in% c('sbm','unconstr')) newdata = design
  yhat1 = plogis(newdata%*%coef(comp1))
  yhat2 = plogis(newdata%*%coef(comp2))
  yhat3 = plogis(newdata%*%coef(comp3))
  yhat4 = plogis(newdata%*%coef(comp4))
  yhat5 = plogis(newdata%*%coef(comp5))
  loglik = logLik(comp1)+logLik(comp2)+logLik(comp3)+logLik(comp4)+logLik(comp5)
  return(list('comp1'=comp1, 'comp2'=comp2, 'comp3'=comp3, 'comp4'=comp4, 'comp5'=comp5,
              'yhat1' = yhat1, 'yhat2' = yhat2, 'yhat3' = yhat3, 'yhat4' = yhat4, 'yhat5' = yhat5,
              'loglik' = as.numeric(loglik)))
}

# glm estimator with 6 components, p1 or sbm models:
glm_estimator_6comp = function(ylist, weights, design, model = c('p1', 'sbm', 'unconstr'), family = 'binomial') {
  if (!exists('v')) v = (1+sqrt(1+8*length(ylist[[1]])))/2
  w1 = rep(weights[,1], each = v*(v-1)/2)
  w2 = rep(weights[,2], each = v*(v-1)/2)
  w3 = rep(weights[,3], each = v*(v-1)/2)
  w4 = rep(weights[,4], each = v*(v-1)/2)
  w5 = rep(weights[,5], each = v*(v-1)/2)
  w6 = rep(weights[,6], each = v*(v-1)/2)
  y = do.call(c, ylist)
  K = length(ylist) # number of graphs
  X = do.call( rbind, replicate(K, design, simplify=FALSE) )  
  
  if (model == 'p1') {
    comp1 = glm(y ~ X, weights = w1, family = family)
    comp2 = glm(y ~ X, weights = w2, family = family)
    comp3 = glm(y ~ X, weights = w3, family = family)
    comp4 = glm(y ~ X, weights = w4, family = family)
    comp5 = glm(y ~ X, weights = w5, family = family)
    comp6 = glm(y ~ X, weights = w6, family = family)
  }
  else if (model %in% c('sbm','unconstr')) {
    comp1 = glm(y ~ X-1, weights = w1, family = family)
    comp2 = glm(y ~ X-1, weights = w2, family = family)
    comp3 = glm(y ~ X-1, weights = w3, family = family)
    comp4 = glm(y ~ X-1, weights = w4, family = family)
    comp5 = glm(y ~ X-1, weights = w5, family = family)
    comp6 = glm(y ~ X-1, weights = w6, family = family)
  }   
  
  if (model == 'p1') newdata = cbind(1, design)
  else if (model %in% c('sbm','unconstr')) newdata = design
  yhat1 = plogis(newdata%*%coef(comp1))
  yhat2 = plogis(newdata%*%coef(comp2))
  yhat3 = plogis(newdata%*%coef(comp3))
  yhat4 = plogis(newdata%*%coef(comp4))
  yhat5 = plogis(newdata%*%coef(comp5))
  yhat6 = plogis(newdata%*%coef(comp6))
  loglik = logLik(comp1)+logLik(comp2)+logLik(comp3)+logLik(comp4)+logLik(comp5)+logLik(comp6)
  return(list('comp1'=comp1, 'comp2'=comp2, 'comp3'=comp3, 'comp4'=comp4, 'comp5'=comp5, 'comp6'=comp6,
              'yhat1' = yhat1, 'yhat2' = yhat2, 'yhat3' = yhat3, 'yhat4' = yhat4, 'yhat5' = yhat5, 'yhat6' = yhat6,
              'loglik' = as.numeric(loglik)))
}

# glm estimator with 7 components, p1 or sbm models:
glm_estimator_7comp = function(ylist, weights, design, model = c('p1', 'sbm', 'unconstr'), family = 'binomial') {
  if (!exists('v')) v = (1+sqrt(1+8*length(ylist[[1]])))/2
  w1 = rep(weights[,1], each = v*(v-1)/2)
  w2 = rep(weights[,2], each = v*(v-1)/2)
  w3 = rep(weights[,3], each = v*(v-1)/2)
  w4 = rep(weights[,4], each = v*(v-1)/2)
  w5 = rep(weights[,5], each = v*(v-1)/2)
  w6 = rep(weights[,6], each = v*(v-1)/2)
  w7 = rep(weights[,7], each = v*(v-1)/2)
  y = do.call(c, ylist)
  K = length(ylist) # number of graphs
  X = do.call( rbind, replicate(K, design, simplify=FALSE) )  
  
  if (model == 'p1') {
    comp1 = glm(y ~ X, weights = w1, family = family)
    comp2 = glm(y ~ X, weights = w2, family = family)
    comp3 = glm(y ~ X, weights = w3, family = family)
    comp4 = glm(y ~ X, weights = w4, family = family)
    comp5 = glm(y ~ X, weights = w5, family = family)
    comp6 = glm(y ~ X, weights = w6, family = family)
    comp7 = glm(y ~ X, weights = w7, family = family)
  }
  else if (model %in% c('sbm','unconstr')) {
    comp1 = glm(y ~ X-1, weights = w1, family = family)
    comp2 = glm(y ~ X-1, weights = w2, family = family)
    comp3 = glm(y ~ X-1, weights = w3, family = family)
    comp4 = glm(y ~ X-1, weights = w4, family = family)
    comp5 = glm(y ~ X-1, weights = w5, family = family)
    comp6 = glm(y ~ X-1, weights = w6, family = family)
    comp7 = glm(y ~ X-1, weights = w7, family = family)
  }   
  
  if (model == 'p1') newdata = cbind(1, design)
  else if (model %in% c('sbm','unconstr')) newdata = design
  yhat1 = plogis(newdata%*%coef(comp1))
  yhat2 = plogis(newdata%*%coef(comp2))
  yhat3 = plogis(newdata%*%coef(comp3))
  yhat4 = plogis(newdata%*%coef(comp4))
  yhat5 = plogis(newdata%*%coef(comp5))
  yhat6 = plogis(newdata%*%coef(comp6))
  yhat7 = plogis(newdata%*%coef(comp7))
  loglik = logLik(comp1)+logLik(comp2)+logLik(comp3)+logLik(comp4)+logLik(comp5)+logLik(comp6)+logLik(comp7)
  return(list('comp1'=comp1, 'comp2'=comp2, 'comp3'=comp3, 'comp4'=comp4, 'comp5'=comp5, 'comp6'=comp6, 'comp7'=comp7,
              'yhat1' = yhat1, 'yhat2' = yhat2, 'yhat3' = yhat3, 'yhat4' = yhat4, 'yhat5' = yhat5, 'yhat6' = yhat6, 'yhat7' = yhat7,
              'loglik' = as.numeric(loglik)))
}
