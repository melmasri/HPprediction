#########################################
## Script to run using cross-validation and uncertainty
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
cleaned = network_clean(com, tree, 'full', uncertainty=TRUE)
com = cleaned$Z                         # cleaned binary interaction matrix
tree = cleaned$tree                     # cleaned tree

## subsetting GMPD for up to 2004 and up to 2010 (com10)
com10 = com
year = 2006
aux = which(com10>year, arr.ind=T)
com = 1*(com10>0)
for(i in 1:nrow(aux))
    if(sum(com[aux[i,1],])>-1 & sum(com[,aux[i,2]])>1) com[aux[i,1], aux[i,2]]<-0
print(sprintf("No. of left out interactions between year %d and end of dataset is %d", year, sum(1*(com10>0)) - sum(com>0)))
print(sprintf("accounts for %f%%  of the data", 100*(sum(1*(com10>0)) - sum(com>0))/sum(com10>0)))
com10=1*(com10>0)


## load useful network analysis functions
source('network_analysis.R')

## indexing 5-folds of interactions
folds = cross.validate.fold(com, n= 5, 1)  # a matrix of 3 columns (row, col, group), (row, col) correspond to Z, group to the CV group
tot.gr = length(unique(folds[,'gr']))   # total number of CV groups

## A loop ran over all CV groups
res = mclapply(1:tot.gr ,function(x, folds, Z, tree, slice, model.type){
    ## Analysis for a single fold
    Z.train = Z
    Z.train[folds[which(folds[,'gr']==x),c('row', 'col')]]<-0

    ## ################################################
    ## running the model of interest with uncertainty
    obj = network_est(Z.train, slices=slice, tree=tree,
        model.type=model.type, uncertainty = TRUE)
    ## Extracting mean posteriors
    y = if(is.matrix(obj$param$y)) rowMeans(obj$param$y) else  mean(obj$param$y)
    w = if(is.matrix(obj$param$w)) rowMeans(obj$param$w) else  mean(obj$param$w)
    eta = if(is.null(obj$param$eta)) 0 else mean(obj$param$eta)
    g = if(is.null(obj$param$g)) 0 else mean(obj$param$g) 

    withG = list(w=w, y=y, eta=eta, g=g, g.sample = obj$param$g)

    ## ################################################
    ## running the model of interest without uncertainty
    obj = network_est(Z.train, slices=slice, tree=tree,
        model.type=model.type)
    ## Extracting mean posteriors
    y = if(is.matrix(obj$param$y)) rowMeans(obj$param$y) else  mean(obj$param$y)
    w = if(is.matrix(obj$param$w)) rowMeans(obj$param$w) else  mean(obj$param$w)
    eta = if(is.null(obj$param$eta)) 0 else mean(obj$param$eta)

    withOutG = list(w=w, y=y, eta=eta)
    ## ################################################
    ## return
    list(withG = withG, withOutG = withOutG)
    
},folds=folds,Z = com, tree=tree, model.type=MODEL, slice = SLICE,
    mc.preschedule = TRUE, mc.cores = min(tot.gr, NO.CORES))

## Averaging mean posterior estimates

## Model with uncertainty
W = rowMeans(sapply(res, function(r) r[['withG']][['w']]))
Y = rowMeans(sapply(res, function(r) r[['withG']][['y']]))
Eta = mean(sapply(res, function(r) r[['withG']][['eta']]))
G = mean(sapply(res, function(r) r[['withG']][['g']]))
G.sample = sapply(res, function(r)r[['withG']][['g.sample']])

if(grepl('(full|dist)', MODEL)){
    distance = 1/cophenetic(rescale(tree, 'EB', Eta))
    diag(distance)<-0
    distance = distance %*% com
    if(grepl('dist', MODEL)) distance[distance==0]<-Inf else 
    distance[distance==0]<-1
}else distance = 1

## Probability matrix with G
P =  1 - exp(-outer(Y, W)*distance)
Pg = G*P/(1-P + G*P)
Pg[com>0]<-P[com>0]

## Model with uncertainty
W = rowMeans(sapply(res, function(r) r[['withOutG']][['w']]))
Y = rowMeans(sapply(res, function(r) r[['withOutG']][['y']]))
Eta = mean(sapply(res, function(r) r[['withOutG']][['eta']]))

if(grepl('(full|dist)', MODEL)){
    distance = 1/cophenetic(rescale(tree, 'EB', Eta))
    diag(distance)<-0
    distance = distance %*% com
    if(grepl('dist', MODEL)) distance[distance==0]<-Inf else 
    distance[distance==0]<-1
}else distance = 1
## Probability matrix with G
Pnog =  1 - exp(-outer(Y, W)*distance)


## left ordering based on com10
indices = lof(com10, indices = TRUE)
com10 = com10[, indices]
com = com[, indices]
Pg = Pg[, indices]
Pnog = Pnog[, indices]

## Histogram of G
pdf(paste0(subDir, 'hist_g.pdf'),height=4)
hist(G.sample,freq=F,col='ivory4', xlab='Posterior estimate of g', main='')
abline(v=quantile(G.sample,probs = c(0.05, 0.95)), col="red", lty=2, lwd=2)
dev.off()
## Printing posterior mean and empirical quantiles
write.csv(cbind(G = mean(G.sample), t(quantile(G.sample,probs = c(0.05, 0.95)))),
              file= paste0(subDir, 'G-posterior-quantile.csv'))
print('Posterior estimate and empirical quantiles are')
print(cbind(G = mean(G.sample), t(quantile(G.sample,probs = c(0.05, 0.95)))))


## printing input Z
pdf(paste0(subDir, 'Z_input_2004.pdf'))
plot_Z(com)
dev.off()

## printing input Z
pdf(paste0(subDir, 'Z_input_2010.pdf'))
plot_Z(com10)
dev.off()

## printing input tree
pdf(paste0(subDir, 'tree_input.pdf'))
plot(tree, show.tip.label=FALSE)
dev.off()

name = 'GMPD'
## ROC curves and posterior interaction matrices
rocNoG = rocCurves(com10, com,Pnog, all=TRUE, plot=FALSE, bins=400)
pdf(paste0(subDir, 'without_gZ-2010_',name,'.pdf'))
plot_Z(1*(Pnog> rocNoG$threshold))
dev.off()

rocG= rocCurves(com10, com, Pg, all=TRUE, plot=FALSE, bins=400)
pdf(paste0(subDir,'with_gZ-2010_',name,'.pdf'))
plot_Z(1*(Pg> rocG$threshold))
dev.off()

## Printing AUC and prediction results
tb = data.frame(rbind(ana.table(com10, com, Pg,   roc = rocG), 
                 ana.table(com10, com, Pnog, roc = rocNoG)))
rownames(tb)=c('Model with $g$', 'Model without $g$')
tb = tb[,c('auc','pred.held.out.ones','pred.tot.ones')]
print(tb)
write.csv(tb, file=paste0(subDir, 'AUC-PRED.tex'))

## Printing ROC curves
pdf(paste0(subDir, 'ROC-g_', MODEL,'.pdf'))
names = paste(MODEL, c('with g', 'without g'))
plot(rocG$roc$FPR, rocG$roc$TPR, xlab='1-specificity', ylab = 'sensitivity',type ='b', col='black', main = 'ROC Curve', xlim = c(0,1), ylim = c(0,1), lty=1, lwd=3, pch=5)
lines(rocNoG$roc$FPR, rocNoG$roc$TPR,type='b', col='ivory4',lty=2, lwd=3,pch=2)
abline(a = 0, b=1,col='black',lty=2, lwd=2)
legend('bottomright', legend = names, col=c('black', 'ivory4'), lty=c(1,2),lwd=3, pch=c(5,2))
dev.off()


## Degree Distribution with and without G
ZpostG = 1*(Pg>rocG$threshold)
ZpostNoG = 1*(Pnog>rocNoG$threshold)
pdf(paste0(subDir,'degree_dist.pdf'))
par(mfrow = c(2,2))
plot_degree(com10, ZpostG, type='hosts',host.col = 'red')
plot_degree(com10, ZpostG, type='parasites', parasite.col = 'blue')
plot_degree(com10, ZpostNoG, type='hosts',host.col = 'red')
plot_degree(com10, ZpostNoG, type='parasites', parasite.col = 'blue')
dev.off()

    
## Printing top M interactions
## top m
m = 4000
## Full with G
ord.pg = order(Pg, decreasing=TRUE)
## Full without G
ord.p = order(Pnog, decreasing=TRUE)

topm = data.frame(t(sapply(1:m, function(r) c(actual = sum(com10[ord.pg[1:r]]),
    withG= sum(com10[ord.pg[1:r]]*ZpostG[ord.pg[1:r]]),
    withOutG= sum(com10[ord.p[1:r]]*ZpostNoG[ord.p[1:r]])))))

pdf(paste0(subDir, 'TopM_.pdf'))
plot(x=1:m,y = topm[,'withG'], xlab='Number of validated pairwise interactions', ylab = 'Number of recovered pairwise interactions',  col='red',lty=1, type='l', lwd=2)
lines(x=1:m,y = topm[,'withOutG'], lty=3, type='l', lwd=2, col='red')
lines(x=1:m, y=1:m, lty=5, lwd=1, col='black')
gnames = c(paste(MODEL, c('with uncertainty', '')), 'x=y')
legend('bottomright', legend = gnames,col=c('red','red','black'),lty=c(1,3,5),lwd=c(2,1))
dev.off()


## Saving workspace
if(SAVE_PARAM)
    save.image(file = paste0(subDir, SAVE_FILE))

##################################################
##################################################
