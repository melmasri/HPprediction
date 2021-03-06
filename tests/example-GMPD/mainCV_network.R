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

## loading mammal supertree included in Fritz et al. (2009)
source('data/download_tree.R')  # see variable 'tree'

## loading GMPD
if(exists("PATH.TO.FILE") && !is.null(PATH.TO.FILE)){
    if(grepl('.rds', PATH.TO.FILE, ignore.case = TRUE))
        com <- readRDS(PATH.TO.FILE) else 
    load(PATH.TO.FILE)
}else{
    source('data/load_GMPD.R')           # see matrix 'com'    
}


## aux = which(colSums(1*(com>0))==1)
## com = com[, -aux]
## com = com[-which(rowSums(1*(com>0))==0), ]
cat("dim com:", dim(com), 'no.interactions:', sum(1*(com>0)), '\n')

## preparing tree and com
cleaned = network_clean(com, tree, 'full')
com = cleaned$Z                         # cleaned binary interaction matrix
tree = cleaned$tree                     # cleaned tree

## indexing 5-folds of interactions
## warning('CV is run using minimum of 2 interactions per column! to change this ')
folds = cross.validate.fold(com, n= 5, CV.MIN.PER.COL)  # a matrix of 3 columns (row, col, group), (row, col) correspond to Z, group to the CV group
tot.gr = length(unique(folds[,'gr']))   # total number of CV groups

## A loop ran over all CV groups
res = mclapply(1:tot.gr ,function(x, folds, Z, tree, slice, model.type, ALPHA.ROWS, ALPHA.COLS){
    ## Analysis for a single fold
    Z.train = Z
    Z.train[folds[which(folds[,'gr']==x),c('row', 'col')]]<-0

    ## running the model of interest
    obj = network_est(Z.train, slices=slice, tree=tree, model.type=model.type,
                      a_y = ALPHA.ROWS, a_w = ALPHA.COLS)

    P = sample_parameter(obj$param, model.type, Z.train, tree)
    Eta = if(is.null(obj$param$eta)) 0 else mean(obj$param$eta)
    
    ## order the rows in Z.test as in Z.train
    roc = rocCurves(Z, Z.train, P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, Z.train, P, roc,  plot=FALSE)
    roc.all = rocCurves(Z, Z.train, P=P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, Z.train, P, roc.all, plot=FALSE)
    
    list(param=list(P=P, Eta = Eta), tb = tb,
         tb.all = tb.all, FPR.all = roc.all$roc$FPR,
         TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
    
},folds=folds,Z = com, tree=tree, model.type=MODEL, slice = SLICE,
    ALPHA.COLS= ALPHA.COLS, ALPHA.ROWS = ALPHA.ROWS,
    mc.preschedule = TRUE, mc.cores = min(tot.gr, NO.CORES))

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
P = matrix(rowMeans(sapply(res, function(r) r$param$P)),
    nrow = nrow(com), ncol = ncol(com))

## left ordering of interaction and probability matrix
indices = lof(com, indices = TRUE)
com = com[, indices]
P = P[, indices]
rownames(P)<-rownames(com)
colnames(P)<-colnames(com)
## print topPairs
topPairs(P,1*(com>0),topX=50)

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

## printing output tree
if(grepl('(full|dist)', MODEL)){
    Eta = mean(sapply(res, function(r) r$param$Eta))
    print(paste('Eta is', Eta))
    pdf(paste0(subDir, 'tree_', MODEL,'.pdf'))
    plot(rescale(tree, 'EB', Eta), show.tip.label=FALSE)
    dev.off()
}
## Saving workspace
if(SAVE_PARAM)
    save.image(file = paste0(subDir, SAVE_FILE))

##################################################
##################################################

