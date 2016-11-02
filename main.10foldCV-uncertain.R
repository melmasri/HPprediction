#########################################
## Script to introduce uncertainty

## Global Variable
SAVE_PARAM = TRUE
# DATAFILENAME = 'comEID-subset.RData'
## DATAFILENAME = 'comGMPD-year.RData'
## DATAFILENAME = 'comGMPD-year.single.RData'
## SUBSET = FALSE
dataset = if(grepl('GMP', DATAFILENAME)) 'gmp' else 'eid'
print(DATAFILENAME)
## source('library.R')
## source('gen.R')
load(DATAFILENAME)
library(parallel)

#######################
## subsetting
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
    aux = rowSums(1*(com>0))
    com = com[aux>1,]
    phy_dist = phy_dist[aux>1,]
    phy_dist = phy_dist[,aux>1]
    com = lof(com)
}
slice = max(ceiling(14000/ncol(com)),5)
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
com_pa = 1*(com>0)
    
if(grepl('GMP', DATAFILENAME)){
    com10 = com
    year = 2004
    aux = which(com10>year, arr.ind=T)
    com = 1*(com10>0)
    for(i in 1:nrow(aux))
        if(sum(com[aux[i,1],])>1 & sum(com[,aux[i,2]])>1)
            com[aux[i,1], aux[i,2]]<-0
    print(sprintf("No. of left out interactions between year %d and end of dataset is %d", year, sum(1*(com10>0)) - sum(com>0)))
    print(sprintf("accounts for %f%%  of the data", 100*(sum(1*(com10>0)) - sum(com>0))/sum(com10>0)))
    com10=1*(com10>0)
}

if(grepl('EID', DATAFILENAME)){
    com = 1*(com>0)
    com10 = 1*(com>0)
    com = cross.validate.set(com10, 0.1)
}

##  hyper = list(parasite =c(0.32, 1), host = c(0.94,1), eta = c(0.012)) 
pdf('all.pdf')

## without G
### It phylogeny does not improve subset 
if(SUBSET){
    param = gibbs_one(com,slice=slice , hyper=hyper,updateHyper=FALSE, AdaptiveMC = TRUE)
}else
    param = gibbs_one(com,slice=slice,dist=phy_dist, eta=1.5, hyper=hyper,updateHyper=FALSE, AdaptiveMC = TRUE)

ana.plot(param)
paramMu = getMean(param)
if(SUBSET){
    P = 1-exp(-outer(paramMu$y, paramMu$w))
}else{
    P = 1-exp(-outer(paramMu$y, paramMu$w)*((phy_dist^paramMu$eta)%*% com))
}
par(mfrow=c(1,1))
roc = rocCurves(Z =1*(com10>0), Z_cross = com, P=P, plot=TRUE, all=FALSE, bins=400)
ana.table(com10, com, roc, plot=TRUE)

roc.all = rocCurves(Z =1*(com10>0), Z_cross = com, P=P, plot=TRUE, all=TRUE, bins=400)
ana.table(com10, com, roc.all, plot=TRUE)

## Adapting hyper parameters
hyperNoG<-hyper
hyperNoG[['etaSamplingSD']]<-param$sd$eta
hyperNoG[['hostSamplingSD']]<-param$sd$y
hyperNoG[['parasiteSamplingSD']]<-param$sd$w

hyperNoG[['hostStart']]<- paramMu$y
hyperNoG[['parasiteStart']]<- paramMu$w
hyperNoG[['etaStart']]<-paramMu$eta

##################################################
### with G
if(SUBSET){
    paramG = gibbs_one(com,slice=slice , hyper=hyper,updateHyper=FALSE, AdaptiveMC = TRUE, uncertain = TRUE)
}else
    paramG = gibbs_one(com,slice=slice,dist=phy_dist, eta=1, hyper=hyper,updateHyper=FALSE, AdaptiveMC = TRUE, uncertain=TRUE)

ana.plot(paramG)

paramMuG = getMean(paramG)
if(SUBSET){
    PG = 1-exp(-outer(paramMu$y, paramMu$w))
}else{
    PG = 1-exp(-outer(paramMu$y, paramMu$w)*((phy_dist^paramMu$eta)%*% com))
}

PG1 = paramMuG$g*PG/(1-PG  + paramMuG$g*PG)
PG1[com==1]<-PG[com==1]

rocG = rocCurves(Z =1*(com10>0), Z_cross = com, P=PG1, plot=TRUE, all=FALSE, bins=400)
cbind(Model=c('with G', 'without G'), rbind( round(ana.table(com10, com, rocG), 4), round(ana.table(com10, com, roc), 4)))

rocG.all = rocCurves(Z =1*(com10>0), Z_cross = com, P=PG1, plot=TRUE, all=TRUE, bins=400)
cbind(Model=c('with G', 'without G'), rbind( round(ana.table(com10, com, rocG.all), 4), round(ana.table(com10, com, roc.all), 4)))

par(mfrow=c(1,1))
plot(cbind(rocG$roc$FPR, rocG$roc$TPR), type='b', col='red', xlab='1-specificity', ylab = 'sensitivity', main = 'ROC Curve', xlim = c(0,1), ylim = c(0,1))
lines(cbind(roc$roc$FPR, roc$roc$TPR), type='b', col='blue')

par(mfrow=c(1,1))
plot(cbind(rocG.all$roc$FPR, rocG.all$roc$TPR), type='b', col='red', xlab='1-specificity', ylab = 'sensitivity', main = 'ROC Curve', xlim = c(0,1), ylim = c(0,1))
lines(cbind(roc.all$roc$FPR, roc.all$roc$TPR), type='b', col='blue')

dev.off()

hyperG<-hyper
hyperG[['etaSamplingSD']]<-paramG$sd$eta
hyperG[['hostSamplingSD']]<-paramG$sd$y
hyperG[['parasiteSamplingSD']]<-paramG$sd$w

hyperG[['hostStart']]<- paramMuG$y
hyperG[['parasiteStart']]<- paramMuG$w
hyperG[['etaStart']]<-paramMuG$eta


##  FOLD CV for uncertainty
pairs = cross.validate.fold(com, n=5)
tot.gr = length(unique(pairs[,'gr']))

res = mclapply(1:tot.gr ,function(x, pairs, Z, dist,hyperNoG, hyperG, SUBSET){
    source('../library.R', local=TRUE)
    source('../gen.R', local=TRUE)

    slice = max(5,ceiling(6000/ncol(Z)))
    com_paCross = Z
    com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
    ## with G
    if(SUBSET){
        param = gibbs_one(com_paCross,slice=slice,hyper=hyperG,updateHyper=FALSE, AdaptiveMC = TRUE, uncertain = TRUE)
    }else
        param = gibbs_one(com_paCross,slice=slice,dist=dist, eta=1, hyper=hyperG,updateHyper=FALSE, AdaptiveMC = TRUE, uncertain=TRUE)
    
    aux  = getMean(param)
    if(SUBSET){
        P1 = 1-exp(-outer(aux$y, aux$w))
    }else{
        P1 = 1-exp(-outer(aux$y, aux$w)*((dist^aux$eta)%*% com_paCross))
    }
    
    P = aux$g*P1/(1-P1  + aux$g*P1)
    P[com_paCross==1]<-P1[com_paCross==1]
    roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
    roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)

    withG = list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)

    ## without G
    if(SUBSET){
        param = gibbs_one(com_paCross,slice=slice , hyper=hyperNoG,updateHyper=FALSE, AdaptiveMC = TRUE, uncertain = FALSE)
    }else
        param = gibbs_one(com_paCross,slice=slice,dist=dist, eta=1, hyper=hyperNoG,updateHyper=FALSE, AdaptiveMC = TRUE, uncertain=FALSE)
    
    aux  = getMean(param)
    if(SUBSET){
        P = 1-exp(-outer(aux$y, aux$w))
    }else{
        P = 1-exp(-outer(aux$y, aux$w)*((dist^aux$eta)%*% com_paCross))
    }
    
    roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
    roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
    
    withOutG = list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)

    ## Dist only
    DistOnlyWithG<-NULL
    DistOnlyWithOutG<-NULL
    if(!SUBSET){
        ## without G
        param = gibbs_one(com_paCross,slice=slice,dist=dist, eta=1, hyper=hyperNoG,updateHyper=FALSE, AdaptiveMC = TRUE, uncertain=FALSE, distOnly=TRUE)
        aux  = getMean(param)
        P = 1-  exp(-(dist^aux$eta) %*% com_paCross)
        roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
        tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
        roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
        tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
        
        DistOnlyWithOutG = list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
        ## with G
        param = gibbs_one(com_paCross,slice=slice,dist=dist, eta=1, hyper=hyperG,updateHyper=FALSE, AdaptiveMC = TRUE, uncertain=TRUE, distOnly=TRUE)
        aux  = getMean(param)
        P1 = 1-  exp(-(dist^aux$eta) %*% com_paCross)
        P = aux$g*P1/(1-P1  + aux$g*P1)
        P[com_paCross==1]<-P1[com_paCross==1]
        
        roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
        tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
        roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
        tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
        
        DistOnlyWithG = list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
    }
    
    ## Affinity only
    affinWithOutG<-NULL
    affinWithG<-NULL
    if(!SUBSET){
        ## with G
        param = gibbs_one(com_paCross,slice=slice,hyper=hyperG,updateHyper=FALSE, AdaptiveMC = TRUE, uncertain=TRUE)
        aux  = getMean(param)
        P1 = 1-exp(-outer(aux$y, aux$w))
        P = aux$g*P1/(1-P1  + aux$g*P1)
        roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
        tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
        roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
        tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
        affinWithG = list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)

        ## withOut G
        param = gibbs_one(com_paCross,slice=slice,hyper=hyperNoG,updateHyper=FALSE, AdaptiveMC = TRUE, uncertain=FALSE)
        aux  = getMean(param)
        P = 1-exp(-outer(aux$y, aux$w))
        roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
        tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
        roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
        tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
        affinWithOutG = list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
    }
    list(withG=withG, withOutG = withOutG, affinWithOutG= affinWithOutG, affinWithG = affinWithG,DistOnlyWithOutG = DistOnlyWithOutG, DistOnlyWithG = DistOnlyWithG )
},pairs=pairs,Z = com,dist=phy_dist,hyperNoG=hyperNoG, hyperG=hyperG, SUBSET = SUBSET,mc.preschedule = TRUE, mc.cores = tot.gr) 

if(SAVE_PARAM)
    save.image(file = 'param.RData')

## ##################################################
## ##################################################
