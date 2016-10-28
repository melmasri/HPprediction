

com_pa = 1*(com>0)
r = which.max(rowSums(com_pa));r
c = which.max(colSums(com_pa));c
burn = param_phy$burn_in - 20000:1
dim(param_phy$y)
burn = burn[burn>0]
range(burn)
aux = getMean(param_phy)

par(mfrow=c(3,1))
plot(param_phy$eta)
abline(h=aux$eta, col='red')
plot(param_phy$w[c,])
abline(h=aux$w[c], col='red')
plot(param_phy$y[r,])
abline(h=aux$y[r], col='red')

par(mfrow=c(3,1))
plot(param_phy$eta)
abline(h=aux$eta, col='red')
plot(param_phy$w[c,])
abline(h=aux$w[c], col='red')
plot(param_phy$y[r,])
abline(h=aux$y[r], col='red')

par(mfrow=c(2,1))
plot(param_phy$hh[1,], ylab='hosts')
plot(param_phy$hh[3,],ylab='parasites')


range(exp(-aux$eta*dd))
range(as.vector(P))
summary(as.vector(P))

plot(param_phy$y[1,], ylim=c(0, 20))
for(i in 1:10)
    lines(param_phy$y[i, ], col=i)


### nodist
aux = getMean(param_phy)
P = 1-  exp(-outer(aux$y, aux$w))
P = 1-exp(-outer(aux$y, aux$w)*((phy_dist^aux$eta)%*% com))
roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=TRUE, bins=400, all=FALSE)
tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
tb.all

 tb
    auc     thresh tot.inter hold.out      pred  pred.all
1 87.85 0.02506266      3730      620 0.7403226 0.7450402
> tb.all
   auc     thresh tot.inter hold.out      pred  pred.all
1 88.9 0.02005013      3730      620 0.7677419 0.7806971

tb
    auc     thresh tot.inter hold.out      pred  pred.all
1 94.28 0.01503759      3730      620 0.8725806 0.8680965
> tb.all
    auc     thresh tot.inter hold.out      pred  pred.all
1 93.78 0.01503759      3730      620 0.8725806 0.8680965
> 
