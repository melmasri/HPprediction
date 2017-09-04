#########################################
## Script to run using cross-validation 
rm(list= ls())
## General variables
## please specifiy the following parameters
SAVE_PARAM = TRUE
SAVE_FILE = 'param.RData'
TYPE = 'full'                           # full, distance or affinity
SLICE = 100                             # no of iterations
NO.CORES = 2                            # maximum cores to use

## Loading required packages
library(ape)
library(geiger)
library(fulltext)
library(parallel)

## loading mammal supertree included in Fritz et al. (2009)
source('example-GMPD/download_tree.R')  # see variable 'tree'

## loading GMPD
source('example-GMPD/load_GMPD.R')      # see matrix 'com'

## sourcing MCMC script
source('network_MCMC.R')

## preparing tree and com
cleaned = network_clean(com, tree, 'full')
com = cleaned$Z                         # cleaned binary interaction matrix
tree = cleaned$tree                     # cleaned tree

## load useful network analysis functions
source('network_analysis.R')

## indexing 5-folds of interactionsx
folds = cross.validate.fold(com, n= 5)  # a matrix of 3 columns (row, col, group), (row, col) correspond to Z, group to the CV group
tot.gr = length(unique(folds[,'gr']))   # total number of CV groups

## A loop ran over all CV groups
res = mclapply(1:tot.gr ,function(x, folds, Z, tree, slice, model.type){
    ## Analysis for a single fold
    Z.train = Z
    Z.train[folds[which(folds[,'gr']==x),c('row', 'col')]]<-0
    eta_sd = 0.005
    a_y = a_w = 0.15
    ## running the model of interest
    obj = network_est(Z.train, slices=slice, tree=tree, model.type=model.type)
    ## Probability matrix
    ## Extracting mean posteriors
    y = if(is.matrix(obj$param$y)) rowMeans(obj$param$y) else  mean(obj$param$y)
    w = if(is.matrix(obj$param$w)) rowMeans(obj$param$w) else  mean(obj$param$w)
    eta = if(is.null(obj$param$eta)) 0 else mean(obj$param$eta)

    if(!is.null(obj$param$eta)){
        distance = 1/cophenetic(rescale(tree, 'EB', eta))
        diag(distance)<-0
        distance = distance %*% Z.train
        if(grepl('dist', model.type)) distance[distance==0]<-Inf else 
        distance[distance==0]<-1
    }else distance = 1

    ## Probability matrix
    if(grepl('dist', model.type))  P = 1-exp(-distance) else
    P  =1 - exp(-outer(y, w)*distance)
    
    ## order the rows in Z.test as in Z.train
    roc = rocCurves(Z, Z.train, P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, Z.train, P, roc,  plot=FALSE)
    roc.all = rocCurves(Z, Z.train, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, Z.train, P, roc.all, plot=FALSE)
    
    list(param=list(y=y, w=w, eta=eta), tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
    
},folds=folds,Z = com, tree=tree, model.type=TYPE, slice = SLICE,
    mc.preschedule = TRUE, mc.cores = min(tot.gr, NO.CORES))

## Some analysis results, AUC, %1 recovered
TB = data.frame(
    m.auc = sapply(res, function(r) r$tb$auc),
    m.pred.held.out.ones = sapply(res,function(r) r$tb$pred.held.out.ones),
    m.thresh = sapply(res, function(r) r$tb$thresh),
    m.hold.out = sapply(res, function(r) r$tb$held.out.ones)
)
TB

print(sprintf('Model: %s, AUC: %f and prect 1 recovered from held out: %f',
              TYPE,mean(TB$m.auc), mean(TB$m.pred.held.out.ones)))

## ROC curve points, can plot as plot(ROCgraph)
ROCgraph = cbind(
    FPR = rowMeans(sapply(res, function(r) r$FPR)),
    TPR = rowMeans(sapply(res, function(r) r$TPR)))

## Constructing the P probability matrix from CV results
if(grepl('(full|aff)', TYPE)){
    W = rowMeans(sapply(res, function(r) r$param$w))
    Y = rowMeans(sapply(res, function(r) r$param$y))
    YW = outer(Y, W)
} else W=Y=YW=1

if(grepl('(full|dist)', TYPE)){
    Eta = mean(sapply(res, function(r) r$param$eta))
    distance = 1/cophenetic(rescale(tree, 'EB', Eta))
    diag(distance)<-0
    distance = distance %*% com
    if(grepl('dist', TYPE)) distance[distance==0]<-Inf else 
    distance[distance==0]<-1
}else distance = 1

## Probability matrix
P = 1-  exp(-YW*distance)

## left ordering of interaction and probability matrix
indices = lof(com, indices = TRUE)
com = com[, indices]
P = P[, indices]

## printing posterior interaction matrix
pdf(paste0('Z_', TYPE, '.pdf'))
plot_Z(1*(P > mean(sapply(res, function(r) r$tb$thres))),
       xlab = 'parasites', ylab = 'hosts')
dev.off()

if(SAVE_PARAM)
    save.image(file = SAVE_FILE)

##################################################
##################################################
