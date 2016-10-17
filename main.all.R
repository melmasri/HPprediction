## Main Script

## Global Variable
SAVE_PARAM = TRUE
DATAFILENAME = 'comGMPD.single.RData'
## DATAFILENAME = 'comGMPD.RData'
#DATAFILENAME = 'comEID-subset.single.RData'
## DATAFILENAME = 'comEID-PS.RData'
#DATAFILENAME = 'comEID-PS.single.RData'
#DATAFILENAME = 'comGMPD-year.single.RData'

print(DATAFILENAME)

source('library.R')
source('gen.R')
load(DATAFILENAME)

## sink('report-test.txt')
#######################
## Parameters for the independent GGP
## set the correct prior.
print('Setting the prior.')
## No uncertain

## if(dataset =='gmp')
##     hyper = list(parasite =c(29.8, 1), host = c(0.24,1), eta = c(0.008)) # comGMPD.single.RData

## if(dataset =='eid')
##     hyper = list(parasite= c(0.5, 1), host =c(0.1, 2), eta = c(0.01))
updateHyper = TRUE
AdaptiveMC = TRUE
hyper = list(parasite= c(0.34, 1), host =c(0.93, 2), eta = c(0.005)) # need to check
slice= ceiling(12000/ncol(com))
slice=8
com=1*(com>0)
param_phy = gibbs_one(Z=1*(com>0),slice=slice,dist=phy_dist, eta=1, hyper=hyper, updateHyper = updateHyper, AdaptiveMC=AdaptiveMC)

param_phy = gibbs_one(Z=1*(com>0),slice=slice, hyper=hyper, updateHyper =updateHyper, AdaptiveMC=AdaptiveMC)
com=1*(com>0)

aux = getMean(param_phy);aux$eta
dd = (phy_dist^aux$eta)%*%com
dd = exp(dd)
range(dd)
## P = 1-  exp(-outer(aux$y, aux$w))
P = 1-  exp(-outer(aux$y, aux$w)*dd)
P = 1-  exp(-outer(aux$y, aux$w)*((phy_dist^aux$eta)%*%com))
P = 1-  exp(-outer(aux$y, aux$w))
range(P)
summary(as.vector(P))
roc = rocCurves(Z=1*(com>0), Z_cross= 1*(com>0), P=P, plot=TRUE, bins=400, all=TRUE);roc$auc
print(roc$auc)

## GMP non-single (pdist)
## 86.25 nodist 
## 90.88 with dist
## 89.83 with dist.right

## EID-PS (pdist)
## 86.27 no dist
## 90.9 wth dist
## 93.86 with new(dist)
## 93.55 with orderd dist

## GMP-Carni (pdist)
## 84.91 no dist
## 82.99 with dist
## 86.12 with old pdist


## if(updateHyper){

##     pdf('hyper-param.pdf')
    
##     par(mfrow=c(2,1))
##     plot(param_phy$hh[1,], type='l', main = '', col='blue', xlab = 'Iteration') 
##     plot(param_phy$hh[3,], type='l', main = '', col='blue', xlab = 'Iteration')

##     dev.off()

##     burn =floor(param_phy$burn_in/1.5):(param_phy$burn_in-1)
##     a = rowMeans(param_phy$hh[,burn])
##     names(a)<-c('host-alpha', 'host-tau', 'parasite-alpha', 'parasite-tau')
##     print(DATAFILENAME)
##     print('Hyper parameters:')
##     print(a)
##     print(param_phy$sd$eta)
    
## }

if(SAVE_PARAM)
    save.image(file = 'param.RData')
##################################################
##################################################

com=1*(com>0)
r = which.max(rowSums(com));r
c = which.max(colSums(com));c

## par(mfrow=c(3,1))
## plot(param_phy$y[r,burn], type='l', main='', col='blue', xlab='Iteration', ylab='Host parameter')
## plot(param_phy$w[c,burn], type='l', main = '', col='blue', xlab = 'Iteration', ylab='Parasite parameter')
## plot(param_phy$eta[burn], type='l', main = '', col='blue', xlab = 'Iteration', ylab='Scaling parameter')



par(mfrow=c(3,1))
plot(param_phy$y[r,], type='l', main='', col='blue', xlab='Iteration', ylab='Host parameter')
plot(param_phy$w[c,], type='l', main = '', col='blue', xlab = 'Iteration', ylab='Parasite parameter')
plot(param_phy$eta, type='l', main = '', col='blue', xlab = 'Iteration', ylab='Scaling parameter')

## x =sapply(1:nrow(param_phy$w), function(r)cor(param_phy$eta, param_phy$w[r,]))
## plot(colSums(com), x)
## plot(colSums(com),rowMeans(param_phy$w) )
## plot(rowSums(com),rowMeans(param_phy$y) )


sink()
q('no')
