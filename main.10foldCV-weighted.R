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
## DATAFILENAME = 'comGMPD.RData'
## ## ## print(DATAFILENAME)
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

## x = 1
## dist= phy_dist
## Z  =com
## dataset= 'eid'

res = mclapply(1:tot.gr ,function(x, pairs, Z, dataset, dist){
    source('../library.R', local=TRUE)
    source('../gen.R', local=TRUE)
    if(dataset =='gmp'){
        a_w = 0.6
        b_w = 1
        a_y = 0.5
        b_y = 1
        a_e = 1.13
        b_e = 1
        eta_sd = 0.022
        w_sd = 0.15
    }
    if(dataset =='eid'){
        a_w = 0.5
        b_w = 1
        a_y = 0.1
        b_y = 2
        a_e = 0.96
        b_e = 1
        eta_sd = 0.003
        w_sd = 0.15
    }

    if(dataset =='gmp')
        Z[Z>2]<-log(1+Z)[Z>2]
    if(dataset =='eid')
        Z=log(Z+1)/2

    com_paCross = Z
    com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
    param_phy = gibbs_one(com_paCross,slice=8 ,dist= dist, eta=1, wMH = TRUE)
    aux = getMean(param_phy)
    P = 1-  exp(-outer(aux$y, aux$w^aux$eta)*((dist^aux$eta)%*% com_paCross))
    Z= 1*(Z>0)
    com_paCross = 1*(com_paCross>0)
    roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
    roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
    list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
},pairs=pairs,Z = com,dataset= dataset, dist=phy_dist, mc.preschedule = TRUE, mc.cores = tot.gr) 

if(SAVE_PARAM)
    save.image(file = 'param.RData')

##################################################
##################################################
