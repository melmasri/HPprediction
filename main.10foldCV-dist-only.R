######################################################
## Script for only Phylogeny
######################################################

## Global Variable
SAVE_PARAM = TRUE
## DATAFILENAME = 'comEID-PS.RData'
## DATAFILENAME = 'comGMPD.RData'
print(DATAFILENAME)
## source('library.R')
## source('gen.R')
load(DATAFILENAME)
library(parallel)

#######################
## Parameters for the independent GGP
## set the correct prior.
## Summary statistics
cnames = colnames(com)
rnames = rownames(com)
colnames(com)<-1:ncol(com)
rownames(com)<-1:nrow(com)
colnames(phy_dist)<-1:ncol(phy_dist)
rownames(phy_dist)<-1:nrow(phy_dist)
com_pa = 1*(com>0)
pairs = cross.validate.fold(com_pa, n=5)
tot.gr = length(unique(pairs[,'gr']))

res = mclapply(1:tot.gr ,function(x, pairs, Z, dist, hyper){
    source('../library.R', local=TRUE)
    source('../gen.R', local=TRUE)

    com_paCross = Z
    com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
    slice = 5000
    
    param = gibbs_one(com_paCross,slice=slice,dist=dist, eta=1, hyper=hyper,updateHyper=FALSE, AdaptiveMC = TRUE, uncertain=FALSE, distOnly=TRUE)
    aux  = getMean(param)
    P = 1-  exp(-(dist^aux$eta) %*% com_paCross)

    roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
    roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
    list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
},pairs=pairs,Z = com_pa, hyper=hyper, dist=phy_dist,mc.preschedule = TRUE, mc.cores = max(5, tot.gr)) 

if(SAVE_PARAM)
    save.image(file = 'param.RData')

##################################################
##################################################
