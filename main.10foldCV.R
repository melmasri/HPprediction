############################################
## Script to run 10fold cross validation

## Global Variable
SAVE_PARAM = TRUE
## DATAFILENAME = 'comEID-PS.single.RData'
## DATAFILENAME = 'comGMPD.single.RData'
print(DATAFILENAME)
## source('library.R')
## source('gen.R')
load(DATAFILENAME)
library(parallel)

#######################
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

res = mclapply(1:tot.gr ,function(x, pairs, Z, dist, hyper){
    source('../library.R', local=TRUE)
    source('../gen.R', local=TRUE)

    slice = max(ceiling(10000/ncol(Z)), 5)
    com_paCross = Z
    com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
    
    param_phy = gibbs_one(com_paCross,slice=slice ,dist= dist,
        eta=1, hyper = hyper, updateHyper=FALSE, AdaptiveMC=FALSE)
    aux = getMean(param_phy)
    
    P = 1-  exp(-outer(aux$y, aux$w)*((dist^aux$eta)%*% com_paCross))
    roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=TRUE, bins=400, all=FALSE)
    tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
    roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=TRUE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
    list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
},pairs=pairs,Z = com_pa, dist=phy_dist,hyper=hyper, mc.preschedule = TRUE, mc.cores = tot.gr) 

if(SAVE_PARAM)
    save.image(file = 'param.RData')

##################################################
