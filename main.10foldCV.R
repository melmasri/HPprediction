# ########################################################################
# Mohamad Elmasri (elmasri.m@gmail.com)
# This is script is supposed to replicate the results of:
#  Caron, F. (2012) Bayesian nonparametric models for bipartite graphics
# ########################################################################
# Project dates:	start January 19, 2015
# 					close ongoing
#########################################
## rm(list= ls())
## Global Variable
SAVE_PARAM = TRUE
## DATAFILENAME = 'comEID-PS.RData'
## DATAFILENAME = 'comGMPD.RData'
## print(DATAFILENAME)
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

res = mclapply(1:tot.gr ,function(x, pairs, Z, dataset, dist, SIMPLERHO, hyper){
    source('../library.R', local=TRUE)
    source('../gen.R', local=TRUE)

    ## if(dataset =='gmp')
    ##     hyper = list(parasite =c(29.8, 1), host = c(0.24,1), eta = c(0.008)) #
        
    ## if(dataset =='eid')
    ##     hyper = list(parasite= c(0.5, 1), host =c(0.1, 2), eta = c(0.01))

    slice = min(ceiling(8000/ncol(Z)), 5)
    com_paCross = Z
    com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
    ## dist  = com_paCross%*%t(com_paCross)
    ## dist  = 1-1/dist
    ## dist[is.infinite(dist)]<-0
    ## dist = dist + 1e-4
    ## diag(dist)<-0

    param_phy = gibbs_one(com_paCross,slice=slice ,dist= dist,
        eta=1, wMH = !SIMPLERHO,yMH=FALSE, yEta=FALSE, hyper = hyper, wEta= !SIMPLERHO, updateHyper=FALSE)
    aux = getMean(param_phy)
    if(SIMPLERHO){
        P = 1-  exp(-outer(aux$y, aux$w)*((dist^aux$eta)%*% com_paCross))
    }else{
        P = 1-  exp(-outer(aux$y, aux$w^aux$eta)*((dist^aux$eta)%*% com_paCross))
    }
    roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
    roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
    list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
},pairs=pairs,Z = com_pa,dataset= dataset, dist=phy_dist,SIMPLERHO=SIMPLERHO,hyper=hyper, mc.preschedule = TRUE, mc.cores = tot.gr) 

if(SAVE_PARAM)
    save.image(file = 'param.RData')

##################################################
##################################################
