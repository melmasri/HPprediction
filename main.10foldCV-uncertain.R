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
library(ape)
library(geiger)
#######################
## subsetting
if(!is.null(SUBSET)){
## EID
    if(dataset=='eid') aux = rownames(com) %in% pan$bionomial[pan$Order=="Rodentia"]
    ##GMP
    if(dataset=='gmp') aux = rownames(com) %in% pan$bionomial[pan$Order=="Carnivora"]

    com =com[aux,]
    aux = colSums(1*(com>0))
    com = com[,aux>1]
    aux = rowSums(1*(com>0))
    com = com[aux>1,]
    com = lof(com)
}

## Loading the tree
tree <- read.tree('../mammals.tre')
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(com)])

# Testing all names in hosts in com exist in tree
if(!all(sapply(rownames(com), function(r) r %in% tree$tip.label))) {
    print('Warning! Not all hosts in com exist in tree. Hosts not found in tree will be removed.')
    com <- com[rownames(com)%in%tree$tip.label,]
}

dd = cophenetic(rescale(tree, 'EB', 0))
host.order <- sapply(rownames(dd), function(r) which(r==rownames(com)))
com = com[host.order,]


## Parameters for 
## set the correct prior.
## Summary statistics
cnames = colnames(com)
rnames = rownames(com)
com = unname(com)
com_pa = 1*(com>0)

if(grepl('GMP', DATAFILENAME)){
    com10 = com
    year = 2004
    aux = which(com10>year, arr.ind=T)
    com = 1*(com10>0)
    for(i in 1:nrow(aux))
        if(sum(com[aux[i,1],])>-1 & sum(com[,aux[i,2]])>1) com[aux[i,1], aux[i,2]]<-0
        print(sprintf("No. of left out interactions between year %d and end of dataset is %d", year, sum(1*(com10>0)) - sum(com>0)))
    print(sprintf("accounts for %f%%  of the data", 100*(sum(1*(com10>0)) - sum(com>0))/sum(com10>0)))
    com10=1*(com10>0)
}



##################################################
##################################################
## Global set-up
beta= 0.5
a_y = a_w =  round(sqrt(-log(1- mean(rowMeans(1*(com>0))))),2); a_y;a_w
eta_sd = 0.005

## Printing to a PDF
pdf('all.pdf')
##################################################
##################################################
## without G

param = ICM_est(Z=com,slice=SLICE,tree=tree, eta = 0,
    burn=0.5, eta_sd = eta_sd, a_w =a_w, a_y = a_y, beta=beta)

ana.plot(param, wait=FALSE)
paramMu = getMean(param)
pdist= 1/cophenetic(rescale(tree, 'EB', paramMu$eta))
diag(pdist)<-0
pdist = pdist%*%com
pdist[pdist==0] <-  1
P = 1- exp(-outer(paramMu$y, paramMu$w)*pdist)

par(mfrow=c(1,1))
roc = rocCurves(Z =1*(com10>0), Z_cross = com, P=P, plot=TRUE, all=FALSE, bins=400)
ana.table(com10, com, roc, plot=TRUE)

roc.all = rocCurves(Z =1*(com10>0), Z_cross = com, P=P, plot=TRUE, all=TRUE, bins=400)
ana.table(com10, com, roc.all, plot=TRUE)

##################################################
### with G

paramG = ICM_est(Z=com,slice=SLICE,tree=tree,eta=0,burn=0.5,
    eta_sd = eta_sd, a_w =a_w, a_y = a_y, beta=beta, g=0)

ana.plot(paramG, wait=FALSE)

paramMuG = getMean(paramG)
pdist= 1/cophenetic(rescale(tree, 'EB', paramMuG$eta))
plot(rescale(tree, 'EB', paramMuG$eta))
paramMuG$eta
diag(pdist)<-0
pdist = pdist%*%com
pdist[pdist==0] <-  1
PG = 1- exp(-outer(paramMuG$y, paramMuG$w)*pdist)

PG1 = paramMuG$g*PG/(1-PG  + paramMuG$g*PG)
PG1[com==1]<-PG[com==1]

par(mfrow=c(1,1))
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

##  FOLD CV for uncertainty
pairs = cross.validate.fold(com, n=5,1)
tot.gr = length(unique(pairs[,'gr']))

res = lapply(1:tot.gr ,function(x, pairs, Z, tree, SLICE){
    source('../library.R', local=TRUE)
    source('../gen.R', local=TRUE)
    library(ape)
    library(geiger)

    beta= 0.5
    a_y = a_w =  round(sqrt(-log(1- mean(rowMeans(1*(Z>0))))),2); a_y;a_w
    eta_sd = 0.005
    com_paCross = Z
    com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0

    ## with G
    param = ICM_est(Z=com_paCross,slice=SLICE,tree=tree,eta=0,burn=0.5,
        eta_sd = eta_sd, a_w =a_w, a_y = a_y, beta=beta, g=0)

    aux  = getMean(param)
    
    pdist= 1/cophenetic(rescale(tree, 'EB', aux$eta))
    diag(pdist)<-0
    pdist = pdist%*%com_paCross
    pdist[pdist==0] <-  1
    P1 = 1- exp(-outer(aux$y, aux$w)*pdist)
    P = aux$g*P1/(1-P1  + aux$g*P1)
    P[com_paCross==1]<-P1[com_paCross==1]

    roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
    roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)

    withG = list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)

    ## without G
    param = ICM_est(Z=com_paCross,slice=SLICE,tree=tree,eta=0,burn=0.5,
        eta_sd = eta_sd, a_w =a_w, a_y = a_y, beta=beta)

    aux  = getMean(param)
    pdist= 1/cophenetic(rescale(tree, 'EB', aux$eta))
    diag(pdist)<-0
    pdist = pdist%*%com_paCross
    pdist[pdist==0] <-  1
    P = 1- exp(-outer(aux$y, aux$w)*pdist)

    roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
    roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)

    withOutG = list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
    
    list(withG=withG, withOutG = withOutG)
},pairs=pairs,Z = com,tree=tree, SLICE=SLICE)

##mc.preschedule = TRUE, mc.cores = min(tot.gr, NO.CORES))

if(SAVE_PARAM)
    save.image(file = 'param.RData')


## ##################################################
