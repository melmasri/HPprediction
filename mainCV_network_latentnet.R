#########################################
## Script to run using cross-validation 
#########################################

## General variables
## please specify the following parameters
## SAVE_PARAM = TRUE                    # should workspace be saved
## SAVE_FILE = 'param.RData'            # name of output R workspace file
## MODEL = 'full'                       # full, distance or affinity
## SLICE = 100                          # no of iterations
## subDir = ''                          # directory to print the results 
## NO.CORES = 2                         # maximum cores to use
## COUNT = TRUE                         # TRUE = count data, FALSE = year of first pub.

## Loading required packages
library(parallel)
library(latentnet)
## loading mammal supertree included in Fritz et al. (2009)
source('example-GMPD/download_tree.R')  # see variable 'tree'

## loading GMPD
source('example-GMPD/load_GMPD.R')      # see matrix 'com'
## aux = which(colSums(1*(com>0))==1)
## com = com[, -aux]
## com = com[-which(rowSums(1*(com>0))==0), ]
dim(com)
## sourcing MCMC script
source('network_MCMC.R')

## preparing tree and com
cleaned = network_clean(com, tree, 'full')
com = cleaned$Z                         # cleaned binary interaction matrix
tree = cleaned$tree                     # cleaned tree

## load useful network analysis functions
source('network_analysis.R')

## indexing 5-folds of interactions
folds = cross.validate.fold(com, n= 5, 2)  # a matrix of 3 columns (row, col, group), (row, col) correspond to Z, group to the CV group
tot.gr = length(unique(folds[,'gr']))   # total number of CV groups

## A loop ran over all CV groups
res = mclapply(1:tot.gr ,function(x, folds, Z, tree, slice, model.type){
    ## Analysis for a single fold
    Z.train = Z
    Z.train[folds[which(folds[,'gr']==x),c('row', 'col')]]<-0
    aux = which(rowSums(Z.train)==0)
    ## running the model of interest
    X = if(length(aux)>0) network(Z.train[-aux,]) else network(Z.train)
    fit<-ergmm(X~ euclidean(d=2)+ rsociality,
               control=ergmm.control(mle.maxit=10,burnin=0,threads=2),verbose=TRUE)
    pred <- predict(fit)

    parasites = which(network.vertex.names(X) %in% colnames(Z))
    hosts = which(network.vertex.names(X) %in% rownames(Z))
    P = matrix(0, nrow(Z), ncol(Z))
    if(length(aux)>0)
        P[-aux,] <- t(pred[parasites, hosts]) else
    P= t(pred[parasites, hosts])
    
    ## order the rows in Z.test as in Z.train
    roc = rocCurves(Z, Z.train, P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, Z.train, P, roc,  plot=FALSE)
    roc.all = rocCurves(Z, Z.train, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, Z.train, P, roc.all, plot=FALSE)
    cat(paste('Group ',x,' done..'))
    list(param=list(fit=fit,  P = P), tb = tb,
         tb.all = tb.all, FPR.all = roc.all$roc$FPR,
         TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
    
},folds=folds,Z = com, tree=tree, model.type=MODEL, slice = SLICE,
    mc.preschedule = TRUE, mc.cores = min(tot.gr, NO.CORES))

if(SAVE_PARAM)
    save.image(file = paste0(subDir, SAVE_FILE))

## Some analysis results, AUC, %1 recovered
TB = data.frame(
    m.auc = sapply(res, function(r) r$tb$auc),
    m.pred.held.out.ones = sapply(res,function(r) r$tb$pred.held.out.ones),
    m.thresh = sapply(res, function(r) r$tb$thresh),
    m.hold.out = sapply(res, function(r) r$tb$held.out.ones)
)
TB
## Printing and writing out average MCMC 
print(sprintf('Model: %s, AUC: %f and percent 1 recovered from held out: %f',
              MODEL,mean(TB$m.auc), mean(TB$m.pred.held.out.ones)))

write.csv(rbind(TB, total = colMeans(TB)), file = paste0(subDir, 'AUC-PRED.csv'))

## ROC curve points, can plot as plot(ROCgraph)
ROCgraph = cbind(
    FPR = rowMeans(sapply(res, function(r) r$FPR)),
    TPR = rowMeans(sapply(res, function(r) r$TPR)))

write.csv(ROCgraph, file = paste0(subDir, 'ROC-xy-points.csv'))

## Constructing the P probability matrix from CV results
aux = rowMeans(sapply(res, function(r) r$param$P))
P = matrix(aux, nrow(com), ncol(com))


## left ordering of interaction and probability matrix
indices = lof(com, indices = TRUE)
com = com[, indices]
P = P[, indices]

## printing posterior interaction matrix
pdf(paste0(subDir, 'Z_', MODEL, '.pdf'))
plot_Z(1*(P > mean(sapply(res, function(r) r$tb$thres))))
dev.off()

## printing input Z
pdf(paste0(subDir, 'Z_input.pdf'))
plot_Z(com)
dev.off()

## printing input tree
pdf(paste0(subDir, 'tree_input.pdf'))
plot(tree, show.tip.label=FALSE)
dev.off()

##################################################
##################################################
