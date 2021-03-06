## Script for the nearest NN

MODEL = 'NN'
COUNT = TRUE
SAVE_FILE = 'param.RData'
sTime = Sys.time()
## Formating the sub-directory name
postdix= 'NOSINGLE-FOLD'
subDir =  paste(toupper('cv'),toupper(MODEL),toupper(postdix) ,format(sTime, "%d-%m-%Hh%M"), sep='-')

## adding sub-extension
subDir = paste0('./', subDir, '/')

dir.create(file.path(subDir))
## report file
reportFile <- paste0(subDir, "report.txt")
## Setting the working directory

## starting the process
sink(reportFile)
print(date())

## rm(list= ls())
## Global Variable
SAVE_PARAM = TRUE

set.seed(23456)
## Loading required packages
library(parallel)

## loading mammal supertree included in Fritz et al. (2009)
source('example-GMPD/download_tree.R')  # see variable 'tree'

## loading GMPD
COUNT=TRUE
source('example-GMPD/load_GMPD.R')      # see matrix 'com'
aux = which(colSums(1*(com>0))==1)
com = com[, -aux]
com = com[-which(rowSums(1*(com>0))==0), ]
dim(com)
## sourcing MCMC script
source('network_MCMC.R')

## preparing tree and com
cleaned = network_clean(com, tree, 'full')
com = cleaned$Z                         # cleaned binary interaction matrix
com = 1*(com>0)
tree = cleaned$tree    

## load useful network analysis functions
source('network_analysis.R')

## indexing 5-folds of interactions
folds = cross.validate.fold(com, n= 5, 2)  # a matrix of 3 columns (row, col, group), (row, col) correspond to Z, group to the CV group
tot.gr = length(unique(folds[,'gr']))   # total number of CV groups

## for the optimal nnk.

get.P.nnk<-function(X, nnk){
    P = X * 0
    dist  = X%*%t(X)
    diag(dist)<--10
    nn = lapply(1:nrow(dist),function(r){
        r = dist[r,]
        shared.parasites = sort(r, decreasing = TRUE)[nnk]
        nei= which(r>=shared.parasites)
    })
    for(j in 1:ncol(X))
        for(i in 1:nrow(X))
            if(X[i,j]==0){
                P[i,j]<-mean(X[nn[[i]],j])
            }else{
                ## removing self interaction in determining distance to hosts
                z = X[i,]
                z[j]=0
                nei = tcrossprod(z, X)
                nei[i]<--10
                shar.para = sort(nei, decreasing = TRUE)[nnk]
                P[i,j]<-mean(X[which(nei>=shar.para), j])
            }
    P
}

search.for.nnk<-function(X, Xorg = NULL){
    min.nnk = 2
    max.nnk = nrow(X) - 1
    qartP = floor((max.nnk - min.nnk)/4)
    midP = qartP*2
    qart2P = qartP*3
    pointlist = c(min.nnk, qartP, midP, qart2P, max.nnk)
    auc= c(0,0,0,0,0)
    repeat{
        q = floor((max.nnk - min.nnk)/4)
        pointlist = c(min.nnk, min.nnk + q,  min.nnk + 2*q, min.nnk + 3*q, max.nnk)
        auc1= sapply(pointlist[which(auc==0)], function(nnk){
            P = get.P.nnk(X, nnk)
            if(is.null(Xorg)){
                roc = rocCurves(X, X, P, plot=FALSE, bins=200, all=TRUE)
            }else{
                roc = rocCurves(Xorg, X, P, plot=FALSE, bins=200, all=FALSE)
            }
            roc$auc
        })
        auc[which(auc==0)]<-auc1
        print(rbind(grid = pointlist, AUC = auc))
        a =  which.max(auc)
        min.nnk = pointlist[max(1,a-1)]
        max.nnk = pointlist[min(5, a+1)]
        auc[1] = auc[max(1, a-1)]
        auc[5] = auc[min(5,a+1)]
        auc[2]=auc[3]=auc[4]=0
        if (abs(max.nnk  - min.nnk)<4) break
    }    
    nnk = pointlist[a]
    print(sprintf('nnk is: %d', nnk))
    return(nnk)
}

search.for.nnk.folds<-function(X){
    fol = cross.validate.fold(X, n= 5, 0)  # a matrix of 3 columns (row, col, group), (row, col) correspond to Z, group to the CV group
    tot.gr = length(unique(fol[,'gr']))   # total number of CV groups
    nn = sapply(1:tot.gr, function(x){
        X.train = X
        X.train[fol[which(fol[,'gr']==x),c('row', 'col')]]<-0
        search.for.nnk(X.train, X)
    })
    mean(nn)

}

## nn.k = search.for.nnk(com)
## cat('nn.k:', nn.k, '\n')
res = lapply(1:tot.gr ,function(x, pairs, Z){
    Z.train = Z
    Z.train[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
    nn.k = search.for.nnk(Z.train)
    ##    nn.k = search.for.nnk.folds(Z.train)
    P = get.P.nnk(Z.train, nn.k)
    roc = rocCurves(Z, Z.train, P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(Z, Z.train, P,roc, plot=FALSE)
    roc.all = rocCurves(Z, Z.train, P, plot=FALSE, bins=400, all=TRUE)
    tb.all  = ana.table(Z, Z.train, P,roc.all, plot=FALSE)
    list(nnk=nn.k, P=P,tb = tb, tb.all = tb.all,
         FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR,
         FPR = roc$roc$FPR, TPR=roc$roc$TPR)
},pairs=folds,Z = com)        


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
aux = rowMeans(sapply(res, function(r) r$P))
P = matrix(aux, nrow(com), ncol(com))
rownames(P)<-rownames(com)
colnames(P)<-colnames(com)

## left ordering of interaction and probability matrix
indices = lof(com, indices = TRUE)
com = com[, indices]
P = P[, indices]

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


## Saving workspace
if(SAVE_PARAM)
    save.image(file = paste0(subDir, SAVE_FILE))

##################################################
##################################################
##-----------------------------------------------------------
eTime = Sys.time()
print(sprintf('End time %s.',format(eTime, "%Y-%m-%d %H:%M:%S")))
## Processing time
print('Total processing time')
print(eTime - sTime)

######################
## Closing sink and reverting work directory.
sink()
system(paste("grep '^[^>+;]'", reportFile, ">", paste0(subDir, "report_clean.txt") ))


q('no')

