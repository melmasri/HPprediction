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
pairs = cross.validate.fold(com_pa)
tot.gr = length(unique(pairs[,'gr']))

res = lapply(1:tot.gr ,function(x, pairs, Z, dist){
    source('../library.R', local=TRUE)
    
    com_paCross = Z
    com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
    
    eta = seq(0,3, by=0.1)
    a = sapply(eta, function(t){
        P = 1-  exp(-(dist^t) %*% com_paCross)
        roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
        roc$auc
    })
    eta = eta[which.max(a)]
    P =  P = 1-  exp(-(dist^eta) %*% com_paCross)
    roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
    roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
    list(param=eta, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR, P = P)
},pairs=pairs,Z = com_pa, dist=phy_dist)

if(SAVE_PARAM)
    save.image(file = 'param.RData')

##################################################
##################################################
