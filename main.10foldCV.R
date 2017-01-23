############################################
## Script to run 10fold cross validation

## Global Variable
SAVE_PARAM = TRUE
## TYPE = 'FULL'
## DATAFILENAME = 'comEID-PS.single.RData'
## DATAFILENAME = 'comEID-PS.RData'
## DATAFILENAME = 'comGMPD.single.RData'
## DATAFILENAME = 'comGMPD.RData'
## DATAFILENAME = 'comSim600x200.RData'
print(DATAFILENAME)
print(TYPE)
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
pairs = cross.validate.fold(1*(com>0), n=5)
tot.gr = length(unique(pairs[,'gr']))
com_pa = if(TYPE == 'WEIGHTED') com else 1*(com>0)

res = mclapply(1:tot.gr ,function(x, pairs, Z, dist, hyper, TYPE){
    source('../library.R', local=TRUE)
    source('../gen.R', local=TRUE)

    ## slice = max(ceiling(10000/ncol(Z)), 5)
    slice=1000
    com_paCross = Z
    com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
    
    if(TYPE == 'DISTONLY'){
        param = ICM_est(com_paCross,slice=slice,dist=dist, eta=1,
            hyper=hyper,AdaptiveMC = TRUE, distOnly=TRUE)
        aux  = getMean(param)
        P = 1-  exp(-(dist^aux$eta) %*% com_paCross)
    }
    if(TYPE == 'AFFINITY'){
        param=ICM_est(com_paCross,slice=slice,hyper=hyper,AdaptiveMC = TRUE)
        aux = getMean(param)
        P = 1-exp(-outer(aux$y,aux$w))
    }
    if(TYPE == 'WEIGHTED'){
        if(!all(range(Z)==c(0,1))) stop('Now a binary Z')
        Z=log(Z+1)/2
        com_paCross = Z
        com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
        param = ICM_est(com_paCross,slice=slice ,dist= dist, eta=1,
            hyper=hyper, AdaptiveMC = TRUE)
        aux = getMean(param)
        P = 1-  exp(-outer(aux$y, aux$w)*((dist^aux$eta)%*% com_paCross))
        Z= 1*(Z>0)
        com_paCross = 1*(com_paCross>0)
    }
    if(TYPE == "10FOLD"){
        param = ICM_est(com_paCross,slice=slice ,dist= dist,eta=1,
            hyper = hyper, AdaptiveMC=TRUE,ICM.HORIZ= FALSE)
        aux = getMean(param)
        P = 1-  exp(-outer(aux$y, aux$w)*((dist^aux$eta)%*%com_paCross))
    }

    roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
    roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
    list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
},pairs=pairs,Z = com_pa, dist=phy_dist,hyper=hyper, TYPE=TYPE, mc.preschedule = TRUE, mc.cores = min(tot.gr, 5)) 

if(SAVE_PARAM)
    save.image(file = 'param.RData')

##################################################
