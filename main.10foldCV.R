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
library(ape)
library(geiger)

if(length(grep('com', ls()))==0)
    stop("no object named 'com' in the data file.")

if(length(grep('tree', ls()))==0)
    stop("no object named 'tree' in the data file.")

tree <- read.tree('../mammals.tre')
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(com)])

# Testing all names in hosts in com exist in tree
if(! all(sapply(rownames(com), function(r) r %in% tree$tip.label))) {
    print('Warning! Not all hosts in com exist in tree. Hosts not found in tree will be removed.')
    com <- com[rownames(com)%in%tree$tip.label,]
}

dd = cophenetic(rescale(tree, 'EB', 0))
host.order <- sapply(rownames(dd), function(r) which(r==rownames(com)))
com = com[host.order,]

#######################
## set the correct prior.
## Summary statistics
cnames = colnames(com)
rnames = rownames(com)
com = unname(com)
pairs = cross.validate.fold(1*(com>0), n=5,2)
tot.gr = length(unique(pairs[,'gr']))

if(TYPE == 'WEIGHTED'){
    if(all(range(com)==c(0,1))) stop('command weighted was passed with a binary Z!')
    com_pa = com  
}else com_pa=  1*(com>0)

res = mclapply(1:tot.gr ,function(x, pairs, Z, tree, hyper, TYPE, ICM.HORIZ, slice){
    source('../library.R', local=TRUE)
    source('../gen.R', local=TRUE)
    library(ape)
    library(geiger)

    com_paCross = Z
    com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0

    eta_sd = 0.005
    a_y = a_w = 0.15
    beta = 0.5
    
    if(TYPE == 'DISTONLY'){
        param = ICM_est(Z=com_paCross,slice=slice, tree=tree,eta=0,
            burn=0.5,eta_sd = eta_sd, distOnly= TRUE, beta =beta)
        aux  = getMean(param)
        pdist= 1/cophenetic(rescale(tree, 'EB', aux$eta))
        diag(pdist)<-0
        pdist = pdist%*%com_paCross
        pdist[pdist==0] <-  Inf
        P = 1- exp(-pdist)
    }
    if(TYPE == 'WEIGHTED'){
        if(all(range(Z)==c(0,1))) stop('A binary Z!')
        Z=log(Z+1)/2
        com_paCross = Z
        com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
         param = ICM_est(Z=com_paCross,slice=slice, tree=tree,eta=0,
            burn=0.5,eta_sd = eta_sd, a_w =a_w, a_y= a_y, beta=beta)
        aux = getMean(param)
        pdist= 1/cophenetic(rescale(tree, 'EB', aux$eta))
        diag(pdist)<-0
        pdist = pdist%*%com_paCross
        pdist[pdist==0] <-  1
        P = 1- exp(-outer(aux$y, aux$w)*pdist)
        Z= 1*(Z>0)
        com_paCross = 1*(com_paCross>0)
    }
    if(TYPE == "FOLD"){
        param = ICM_est(Z=com_paCross,slice=slice, tree=tree,eta=0,
            burn=0.5,eta_sd = eta_sd, a_w =a_w, a_y = a_y, beta=beta)
        aux = getMean(param)
        pdist= 1/cophenetic(rescale(tree, 'EB', aux$eta))
        diag(pdist)<-0
        pdist = pdist%*%com_paCross
        pdist[pdist==0] <-  1
        P = 1- exp(-outer(aux$y, aux$w)*pdist)
    }

    roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
    roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
    list(param=aux, tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
},pairs=pairs,Z = com_pa, tree=tree, TYPE=TYPE, slice = SLICE,
    mc.preschedule = TRUE, mc.cores = min(tot.gr, NO.CORES)) 

if(SAVE_PARAM)
    save.image(file = 'param.RData')

##################################################
