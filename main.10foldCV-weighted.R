#########################################
## Script to use counts

## Global Variable
SAVE_PARAM = TRUE
## DATAFILENAME = 'comEID-PS.RData'
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
pairs = cross.validate.fold(com_pa, n=5)
tot.gr = length(unique(pairs[,'gr']))

res = mclapply(1:tot.gr ,function(x, pairs, Z, dist, hyper){
    source('../library.R', local=TRUE)
    source('../gen.R', local=TRUE)
    
    Z=log(Z+1)/2
    ## if(dataset =='gmp')
    ##     Z[Z>2]<-log(1+Z)[Z>2]
    ## if(dataset =='eid')
    ##     Z=log(Z+1)/2
    
    slice = max(5,ceiling(8000/ncol(Z)))
    com_paCross = Z
    com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
    param_phy = gibbs_one(com_paCross,slice=slice ,dist= dist, eta=1, hyper=hyper,updateHyper=FALSE, AdaptiveMC = TRUE)
    
    aux = getMean(param_phy)
    P = 1-  exp(-outer(aux$y, aux$w)*((dist^aux$eta)%*% com_paCross))
    
    Z= 1*(Z>0)
    com_paCross = 1*(com_paCross>0)
    roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
    roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
    list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
},pairs=pairs,Z = com, dist=phy_dist, hyper=hyper, mc.preschedule = TRUE, mc.cores = tot.gr) 

if(SAVE_PARAM)
    save.image(file = 'param.RData')

##################################################
##################################################
