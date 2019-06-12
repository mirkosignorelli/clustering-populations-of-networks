
library(igraph)

load('results/res-simJ.RData')

n.rep = 100

aic = bic.new = matrix(NA, n.rep, ncol = 4)

for (h in 1:n.rep) {
  for (l in 1:n.subcases) {
    v = 20
    K.each = c(10, 30, 60, 100)[l]
    ncomp = 3
    K = ncomp*K.each
    loglvec = loglik.array[h, , l]
    n.param = v*(v-1)/2*seq(1,5)
    # AIC
    temp1 = -2*loglvec + 2*n.param
    aic[h, l] = which(temp1 == min(temp1))
    # BIC
    temp3 = -2*loglvec + log(K)*n.param
    bic.new[h, l] = which(temp3 == min(temp3))
  }
  print(h)
}

apply(aic, 2, table)
apply(bic.new, 2, table)

xlab = 3*c(10, 30, 60, 100)

pdf('res-simJ.pdf', height = 8, width = 16)
par(mfrow = c(2,4))
for (i in 1:4) { 
  title = paste('AIC, K =', xlab[i])
  pmf(aic[,i], title = title, lwd = 3, cex.main = 2.3, cex.axis = 2.3,
           xlab = 'M', xlim = c(1,5))
}
for (i in 1:4) { 
  title = paste('BIC, K =', xlab[i])
  pmf(bic.new[,i], title = title,  lwd = 3, cex.main = 2.3, cex.axis = 2.3,
      xlab = 'M', xlim = c(1,5))
}
dev.off()