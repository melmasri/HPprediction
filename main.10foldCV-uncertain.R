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
#SUBSET=TRUE
#DATAFILENAME = '../comGMPD-year.RData'
#DATAFILENAME = '../comEID-subset.RData'

if(grepl('GMP', DATAFILENAME)) dataset='gmp'
if(grepl('EID', DATAFILENAME)) dataset='eid'

print(DATAFILENAME)
#source('library.R')
#source('gen.R')
load(DATAFILENAME)
library(parallel)

#######################
## subsetting
slice=12
if(SUBSET){
## EID
    if(dataset=='eid') aux = rownames(com) %in% pan$bionomial[pan$Order=="Rodentia"]
##GMP
    if(dataset=='gmp') aux = rownames(com) %in% pan$bionomial[pan$Order=="Carnivora"]
    com =com[aux,]
    phy_dist = phy_dist[aux,]
    phy_dist = phy_dist[,aux]
    aux = colSums(1*(com>0))
    com = com[,aux>1]
    if(dataset=='eid') slice = 40
}
	
diag(phy_dist)<-0


## Parameters for 
## set the correct prior.
## Summary statistics
cnames = colnames(com)
rnames = rownames(com)
colnames(com)<-1:ncol(com)
rownames(com)<-1:nrow(com)
colnames(phy_dist)<-1:ncol(phy_dist)
rownames(phy_dist)<-1:nrow(phy_dist)

if(grepl('GMP', DATAFILENAME)){
    com10 = com
    year = 2005
    aux = which(com10>year, arr.ind=T)
    com = 1*(com10>0)
    for(i in 1:nrow(aux))
        if(sum(com[aux[i,1],])>1 & sum(com[,aux[i,2]])>1)
            com[aux[i,1], aux[i,2]]<-0
    print(sum(1*(com10>0)) - sum(com>0))
    print((sum(1*(com10>0)) - sum(com>0))/sum(com10>0))
    com10=1*(com10>0)
}

if(grepl('EID', DATAFILENAME)){
    com10 = 1*(com>0)
    com = cross.validate.set(com10, 0.1)
}

pairs = cross.validate.fold(com)
tot.gr = length(unique(pairs[,'gr']))

## No uncertain
if(dataset =='gmp')
    hyper = list(parasite= c(1/3, 1), host =c(2, 1), eta = c(0.01))

if(dataset =='eid')
    hyper = list(parasite= c(0.5, 1), host =c(0.1, 2), eta = c(0.01))


pdf('all.pdf')
## without G
paramRegular = gibbs_one(com,slice=slice,dist= phy_dist, eta=1,wMH = TRUE, hyper=hyper)

ana.plot(paramRegular, com)

paramMu = getMean(paramRegular)
if(SIMPLERHO){
    PRegular = 1-exp(-outer(paramMu$y, paramMu$w)*((phy_dist^paramMu$eta)%*% com))
}else{
    PRegular = 1-exp(-outer(paramMu$y, paramMu$w^paramMu$eta)*((phy_dist^paramMu$eta)%*% com))
}

rocRegular = rocCurves(Z =1*(com10>0), Z_cross = com, P=PRegular, plot=TRUE, all=FALSE, bins=400)
ana.table(com10, com, rocRegular, TRUE)

rocRegular.all = rocCurves(Z =1*(com10>0), Z_cross = com, P=PRegular, plot=TRUE, all=TRUE, bins=400)
ana.table(com10, com, rocRegular.all, TRUE)

##################################################
### with G
paramRegularG = gibbs_one(com,slice=slice,dist= phy_dist,eta=1,wMH=TRUE,uncertain =TRUE, hyper = hyper)

ana.plot(paramRegularG, com)

paramMuG = getMean(paramRegularG)
if(SIMPLERHO){
    PRegularG = 1-  exp(-outer(paramMuG$y,paramMuG$w)*((phy_dist^paramMuG$eta)%*% com))
}else{
    PRegularG = 1-  exp(-outer(paramMuG$y,paramMuG$w^paramMuG$eta)*((phy_dist^paramMuG$eta)%*% com))
}
## Method 1
P=PRegularG

rocRegularG = rocCurves(Z =1*(com10>0), Z_cross = com, P=P, plot=TRUE, all=FALSE, bins=400)
ana.table(com10, com, rocRegularG, TRUE)

P = paramMuG$L[1]*PRegularG/(1-PRegularG  + paramMuG$L[1]*PRegularG)
P[com==1]<-PRegularG[com==1]

rocRegularG.all = rocCurves(Z =1*(com10>0), Z_cross = com, P=P, plot=TRUE, all=TRUE, bins=400)
ana.table(com10, com, rocRegularG.all, TRUE)

plot(cbind(rocRegularG.all$roc$FPR, rocRegularG.all$roc$TPR), type='b', col='red')
lines(cbind(rocRegular.all$roc$FPR, rocRegular.all$roc$TPR), type='b', col='blue')
dev.off()
## with uncertain
res = mclapply(1:tot.gr ,function(x, pairs, Z, dist, dataset,s, SIMPLERHO){
    source('../library.R', local=TRUE)
    source('../gen.R', local=TRUE)
    if(dataset =='gmp')
        hyper = list(parasite= c(1/3, 1), host =c(2, 1), eta = c(0.01))

    if(dataset =='eid')
        hyper = list(parasite= c(0.5, 1), host =c(0.1, 2), eta = c(0.01))

    com_paCross = Z
    com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0

    param_phy=gibbs_one(com_paCross,slice=s,dist= dist,eta=1,wMH=!SIMPLERHO,uncertain=TRUE, hyper=hyper,wEta=!SIMPLERHO)
    aux = getMean(param_phy)

    if(SIMPLERHO){
        P1 = 1-  exp(-outer(aux$y, aux$w)*((dist^aux$eta)%*% com_paCross))
    }else{
        P1 = 1-  exp(-outer(aux$y, aux$w^aux$eta)*((dist^aux$eta)%*% com_paCross))
    }
    P = aux$L[1]*P1/(1-P1  + aux$L[1]*P1)
    P[com_paCross==1]<-P1[com_paCross==1]
    roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
    roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
    list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
},pairs=pairs,Z = com,dist=phy_dist, dataset=dataset,s=slice,SIMPLERHO=SIMPLERHO,mc.preschedule = TRUE, mc.cores = tot.gr) 

if(SAVE_PARAM)
    save.image(file = 'param.RData')

##################################################
##################################################
