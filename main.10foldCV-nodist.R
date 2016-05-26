# ########################################################################
# Mohamad Elmasri (elmasri.m@gmail.com)
# This is script is supposed to replicate the results of:
#  Caron, F. (2012) Bayesian nonparametric models for bipartite graphics
# ########################################################################
# Project dates:	start January 19, 2015
# 					close ongoing
#########################################
##rm(list= ls())
## Global Variable
SAVE_PARAM = TRUE
#DATAFILENAME = 'comGMPD.RData'
#DATAFILENAME = '../comGMPD-year.RData'
#DATAFILENAME = '../comEID-subset.RData'

if(grepl('GMP', DATAFILENAME)) dataset='gmp'
if(grepl('EID', DATAFILENAME)) dataset='eid'

print(DATAFILENAME)
#source('library.R')
#source('gen.R')
load(DATAFILENAME)
library(parallel)

## set the correct prior.
## Summary statistics
cnames = colnames(com)
rnames = rownames(com)
colnames(com)<-1:ncol(com)
rownames(com)<-1:nrow(com)
colnames(phy_dist)<-1:ncol(phy_dist)
rownames(phy_dist)<-1:nrow(phy_dist)

pairs = cross.validate.fold(com)
tot.gr = length(unique(pairs[,'gr']))

res = mclapply(1:tot.gr ,function(x, pairs, Z, dataset, hyper){
    source('../library.R', local=TRUE)
    source('../gen.R', local=TRUE)

        
    ## if(dataset =='gmp')
    ##     hyper = list(parasite =c(29.8, 1), host = c(0.24,1), eta = c(0.008)) #

    
    ## if(dataset =='eid')
    ##     hyper = list(parasite= c(0.5, 1), host =c(0.1, 2), eta = c(0.01))

    slice = ceiling(12000/nrow(Z))
    Z=1*(Z>0)
    com_paCross = Z
    com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
    param_phy=gibbs_one(com_paCross,slice=slice,wMH=FALSE, hyper=hyper, wEta=FALSE,AdaptiveMC = FALSE,updateHyper = FALSE, yEta=FALSE,  yMH=FALSE)
    aux = getMean(param_phy)
    P = 1-exp(-outer(aux$y,aux$w))
    roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
    roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
    list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
},pairs=pairs,Z = com, dataset=dataset,hyper=hyper, mc.preschedule = TRUE, mc.cores = tot.gr) 

if(SAVE_PARAM)
    save.image(file = 'param.RData')

##################################################
##################################################
