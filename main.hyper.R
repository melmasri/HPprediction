## Main Script
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
##library scripts

DATAFILENAME = args[1]
## Global Variable
SAVE_PARAM = TRUE

## DATAFILENAME = 'comGMPD.single.RData'
## DATAFILENAME = 'comGMPD.RData'
## DATAFILENAME = 'comEID-subset.single.RData'
## DATAFILENAME = 'comEID-PS.RData'
## DATAFILENAME = 'comEID-PS.single.RData'
## DATAFILENAME = 'comGMPD-year.single.RData'

print(DATAFILENAME)

source('library.R')
source('gen.R')
load(DATAFILENAME)

sink(paste0('hyper-',DATAFILENAME, '.txt'))
#######################
## Parameters for the independent GGP
## set the correct prior.
print('Setting the prior.')
## No uncertain

com=1*(com>0)
updateHyper = TRUE
AdaptiveMC = FALSE

hyper = list(parasite= c(0.34, 1), host =c(0.93, 2), eta = c(0.005)) # need to check
slice= ceiling(8000/ncol(com))

param_phy = gibbs_one(Z=com,slice=slice,dist=phy_dist, eta=1, hyper=hyper, updateHyper =
                          updateHyper, AdaptiveMC=AdaptiveMC)

aux = getMean(param_phy);aux$eta
P = 1-  exp(-outer(aux$y, aux$w)*((phy_dist^aux$eta)%*%com))
pdf(paste0('Plots-', DATAFILENAME, '.pdf'))
roc = rocCurves(Z=1*(com>0), Z_cross= com, P=P, plot=TRUE, bins=400, all=TRUE)
print(roc$auc)

## GMP non-single (pdist)
## 86.25 nodist 
## 90.88 with dist
## 89.83 with dist.right

## EID-PS (pdist)
## 86.27 no dist
## 90.9 wth dist
## 93.86 with new(dist)

## GMP-Carni (pdist)
## 84.91 no dist
## 82.99 with dist
## 86.12 with old pdist

par(mfrow=c(2,1))
plot(param_phy$hh[1,], type='l', main = '', col='blue', xlab = 'Iteration') 
plot(param_phy$hh[3,], type='l', main = '', col='blue', xlab = 'Iteration')


r = which.max(rowSums(com));r
c = which.max(colSums(com));c

par(mfrow=c(3,1))
plot(param_phy$y[r,], type='l', main='', col='blue', xlab='Iteration', ylab='Host parameter')
plot(param_phy$w[c,], type='l', main = '', col='blue', xlab = 'Iteration', ylab='Parasite parameter')
plot(param_phy$eta, type='l', main = '', col='blue', xlab = 'Iteration', ylab='Scaling parameter')

dev.off()

burn =floor(param_phy$burn_in/1.5):(param_phy$burn_in-1)
a = rowMeans(param_phy$hh[,burn])
names(a)<-c('host-alpha', 'host-tau', 'parasite-alpha', 'parasite-tau')
print(DATAFILENAME)
print('Hyper parameters:')
print(a)
print(param_phy$sd$eta)
    


if(SAVE_PARAM)
    save.image(file = paste0('param-',DATAFILENAME))
##################################################
##################################################

sink()
q('no')
