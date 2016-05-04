# ########################################################################
## Mohamad Elmasri (elmasri.m@gmail.com)
## A 100 simulation each to get the average of a 10-fold CV set
## Using the GMPD dataset
# ########################################################################

## rm(list= ls())
## Global Variable
SAVE_PARAM = TRUE
DATAFILENAME = '../comGMPD.RData'
print(DATAFILENAME)
load(DATAFILENAME)
library(parallel)
print(getwd())
#######################
## Parameters for the independent GGP
sim100 = mclapply(1:10, function(i, com, phy_dist){
    source('../library.R', local=TRUE)
    source('../gen.R', local=TRUE)
    cnames = colnames(com)
    rnames = rownames(com)
    colnames(com)<-1:ncol(com)
    rownames(com)<-1:nrow(com)
    colnames(phy_dist)<-1:ncol(phy_dist)
    rownames(phy_dist)<-1:nrow(phy_dist)
    com_pa = 1*(com>0)
    pairs = cross.validate.fold(com_pa, n=4)
    tot.gr = length(unique(pairs[,'gr']))
    dataset='gmp'
    if(dataset =='gmp'){
        a_w = 0.6
        b_w = 1
        a_y = 0.5
        b_y = 1
        a_e = 1.13
        b_e = 1
        eta_sd = 0.007
        w_sd = 0.15
    }
    com = 1*(com>0)
    res = lapply(1:tot.gr ,function(x, pairs, Z, dist){
        com_paCross = Z
        com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
        param_phy = gibbs_one(com_paCross,slice=2 ,dist= dist, eta=1, wMH = TRUE)
        aux = getMean(param_phy)
        P1 = 1-  exp(-outer(aux$y, aux$w^aux$eta)*((dist^aux$eta)%*% com_paCross))
        roc = rocCurves(Z =Z, Z_cross = com_paCross, P=P1, plot=FALSE, all=FALSE, bins=200)
        zz = 1*(P1>roc$threshold)
        predicted = sum(zz[Z==1 & com_paCross==0])/sum(abs(com_paCross - Z)[Z==1])
        tot.inter = sum(Z)
        hold.out = sum(abs(com_paCross - Z)[Z==1])
        list(FPR = roc$roc$FPR,TPR = roc$roc$TPR,
             thres = roc$threshold, auc = roc$auc, hold.out = hold.out, pred = predicted,
             tot.inter = tot.inter, param=aux)
    },pairs=pairs,Z = com_pa, dist=phy_dist)
    FPR = rowMeans(sapply(res, function(r) r$FPR))
    TPR = rowMeans(sapply(res, function(r) r$TPR))
    m.auc= mean(sapply(res, function(r) r$auc))
    m.thresh=mean(sapply(res, function(r) r$thres))
    m.pred=mean(sapply(res,function(r) r$pred))
    m.hold.out=mean(sapply(res, function(r) r$hold.out))
    W = rowMeans(sapply(res, function(r) r$param$w))
    Y = rowMeans(sapply(res, function(r) r$param$y))
    Eta = mean(sapply(res, function(r) r$param$eta))
    TB = data.frame(m.auc, m.thresh, m.pred, m.hold.out)
    list(graph=cbind(FPR, TPR), ana=TB, w=W, y=Y, eta=Eta)
}, com=com, phy_dist= phy_dist, mc.preschedule = TRUE, mc.cores = 20)

if(SAVE_PARAM)
    save.image(file = 'param.RData')

##################################################
##################################################
