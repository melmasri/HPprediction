############################################
## Script to run 10fold cross validation

## Global Variable
SAVE_PARAM = TRUE
## TYPE = 'FULL'
##  DATAFILENAME = 'comEID-PS.single.RData'
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

if(length(grep('com', ls()))==0)
    stop("no object named 'com' in the data file.")

if(length(grep('phy_dist', ls()))==0)
    stop("no object named 'phy_dist' in the data file.")

if(is.null(dim(com)) | is.null(dim(phy_dist)))
    stop("either 'com' or 'phy_dist' are not a matrix type.")

if(!isSymmetric(phy_dist))
    stop("matrix 'phy_dist' is not symmetric.")

if(nrow(com)!= nrow(phy_dist))
    stop("matrix 'phy_dist' doesn't have the same number of rows as 'com'.")

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
if(TYPE == 'WEIGHTED'){
    if(all(range(com)==c(0,1))) stop('command weighted was passed with a binary Z!')
    com_pa = com  
}else com_pa=  1*(com>0)

res = mclapply(1:tot.gr ,function(x, pairs, Z, dist, hyper, TYPE, ICM.HORIZ, slice){
    source('../library.R', local=TRUE)
    source('../gen.R', local=TRUE)

    com_paCross = Z
    com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
    
    if(TYPE == 'DISTONLY'){
        param = ICM_est(com_paCross,slice=slice,dist=dist, eta=1,
            hyper=hyper,AdaptiveMC = TRUE, distOnly=TRUE, ICM.HORIZ = ICM.HORIZ)
        aux  = getMean(param)
        P = 1-  exp(-(dist^aux$eta) %*% com_paCross)
    }
    if(TYPE == 'AFFINITY'){
        param=ICM_est(com_paCross,slice=slice,hyper=hyper,AdaptiveMC = TRUE, ICM.HORIZ = ICM.HORIZ)
        aux = getMean(param)
        P = 1-exp(-outer(aux$y,aux$w))
    }
    if(TYPE == 'WEIGHTED'){
        if(!all(range(Z)==c(0,1))) stop('Not a binary Z!')
        Z=log(Z+1)/2
        com_paCross = Z
        com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
        param = ICM_est(com_paCross,slice=slice ,dist= dist, eta=1,
            hyper=hyper, AdaptiveMC = TRUE,ICM.HORIZ = ICM.HORIZ)
        aux = getMean(param)
        P = 1-  exp(-outer(aux$y, aux$w)*((dist^aux$eta)%*% com_paCross))
        Z= 1*(Z>0)
        com_paCross = 1*(com_paCross>0)
    }
    if(TYPE == "10FOLD"){
        param = ICM_est(com_paCross,slice=slice ,dist= dist,eta=1,
            hyper = hyper, AdaptiveMC=TRUE,ICM.HORIZ= ICM.HORIZ)
        aux = getMean(param)
        P = 1-  exp(-outer(aux$y, aux$w)*((dist^aux$eta)%*%com_paCross))
    }

    roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
    roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
    list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
},pairs=pairs,Z = com_pa, dist=phy_dist,hyper=hyper, TYPE=TYPE,
    ICM.HORIZ = ICM.HORIZ, slice = SLICE,
    mc.preschedule = TRUE, mc.cores = min(tot.gr, NO.CORES)) 

if(SAVE_PARAM)
    save.image(file = 'param.RData')

##################################################
