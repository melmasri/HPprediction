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
if(exists("PATH.TO.FILE") && !is.null(PATH.TO.FILE)){
    if(grepl('.rds', PATH.TO.FILE, ignore.case = TRUE))
        com <- readRDS(PATH.TO.FILE) else 
    load(PATH.TO.FILE)
}else{
    source('example-GMPD/load_GMPD.R')           # see matrix 'com'    
}
## aux = which(colSums(1*(com>0))==1)
## com = com[, -aux]
## com = com[-which(rowSums(1*(com>0))==0), ]

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


## indexing 5-folds of interactions
folds = cross.validate.fold(com, n= 5, CV.MIN.PER.COL)  # a matrix of 3 columns (row, col, group), (row, col) correspond to Z, group to the CV group
tot.gr = length(unique(folds[,'gr']))   # total number of CV groups

## A loop ran over all CV groups
res = mclapply(1:tot.gr ,function(x, folds, Z, tree, slice, model.type, ALPHA.ROWS, ALPHA.COLS){
    ## Analysis for a single fold
    Z.train = Z
    Z.train[folds[which(folds[,'gr']==x),c('row', 'col')]]<-0

    ## ################################################
    ## running the model of interest with uncertainty
    obj = network_est(Z.train, slices=slice, tree=tree,
        model.type=model.type, uncertainty = TRUE, a_y = ALPHA.ROWS, a_w = ALPHA.COLS)
    Eta = if(is.null(obj$param$eta)) 0 else mean(obj$param$eta)
    g = if(is.null(obj$param$g)) 0 else mean(obj$param$g) 
    weights = obj$param$g[sample.int(length(obj$param$g), 1000, replace = TRUE)]
    P = sample_parameter(obj$param, model.type, Z.train, tree, weights  = weights)
    withG = list(P = P, eta=Eta, g=g, g.sample = obj$param$g)

    ## ################################################
    ## running the model of interest without uncertainty
    obj = network_est(Z.train, slices=slice, tree=tree,
        model.type=model.type, a_y = 6, a_w = 0.03)
    ## Extracting mean posteriors
    Pnog = sample_parameter(obj$param, model.type, Z.train, tree)
    Eta = if(is.null(obj$param$eta)) 0 else mean(obj$param$eta)

    withOutG = list(P = Pnog, eta=Eta)
    ## ################################################
    ## return
    list(withG = withG, withOutG = withOutG)
    
},folds=folds,Z = com, tree=tree, model.type=MODEL, slice = SLICE,
    ALPHA.COLS = ALPHA.COLS, ALPHA.ROWS = ALPHA.ROWS,
    mc.preschedule = TRUE, mc.cores = min(tot.gr, NO.CORES))

## Averaging mean posterior estimates

## Model with uncertainty
Pg = matrix(rowMeans(sapply(res, function(r) r[['withG']][['P']])),  nrow = nrow(com), ncol = ncol(com))
G = mean(sapply(res, function(r) r[['withG']][['g']]))
G.sample = sapply(res, function(r)r[['withG']][['g.sample']])

## Model with uncertainty
Pnog = matrix(rowMeans(sapply(res, function(r)  r[['withOutG']][['P']])),  nrow = nrow(com), ncol = ncol(com))

## left ordering based on com10
indices = lof(com10, indices = TRUE)
com10 = com10[, indices]
com = com[, indices]
Pg = Pg[, indices]
Pnog = Pnog[, indices]

## Histogram of G
pdf(paste0(subDir, 'hist_g.pdf'),height=4)
par(mar = c(5,5,1,1)+0.1)
hist(rowMeans(G.sample),freq=F,col='ivory4', xlab='Posterior estimate of g', main='', breaks=20, cex.lab=2, cex.axis = 1.5)
abline(v=quantile(rowMeans(G.sample),probs = c(0.05, 0.95)), col="red", lty=2, lwd=2)
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
names = paste('LS-net: ', MODEL, c('with g', 'without g'))
par(mar = c(5,5,1,1)+0.1)
plot(rocG$roc$FPR, rocG$roc$TPR, xlab='1-specificity', ylab = 'sensitivity',
     type ='l', col='black', main = '', xlim = c(0,1), ylim = c(0,1), lty=1, lwd=3, cex.lab=2, cex.axis=1.5,
      lab=c(x=2,y=2, llen=15))
lines(rocNoG$roc$FPR, rocNoG$roc$TPR,type='l', lty=2, lwd=3)
## abline(a = 0, b=1,col='black',lty=2, lwd=2)
legend('bottomright', legend = names, col=c('black', 'black'), lty=c(1,2),lwd=3, cex=1.5, pt.cex=1.5)
dev.off()

## Printing  obs unk graphs
pdf(paste0(subDir,'with_g_hist_obs_unk_', name,'.pdf'))
colass = rgb(0,0,0,0.3)
colnoass = rgb(0,0,0,0.6)
## colass = rgb(0,0.8,0.8,0.5)
## colnoass = rgb(1,0,0.4,0.4,0.5)
par(mar = c(5,5,1,1)+0.1)
hist(log(Pg[com10>0]),col=colass,main ='', ylab='Density', freq=FALSE, xlab='Log of probability',ylim=c(0, 0.45), xlim=c(-12,0),
     breaks=30, cex.lab=2, lty=1, lwd=4, cex.axis= 1.5)
hist(log(Pg[com10==0]),col=colnoass,freq=FALSE, add=TRUE, lty=2, lwd=4)
legend(x='top', legend=c('Observed associations', 'Unobserved associations'), lwd=4,
       col=c(colass, colnoass),pt.cex=1.5, cex=1.5, lty=c(1,2)) 
dev.off()

pdf(paste0(subDir,'hist_obs_unk_', name,'.pdf'))
colass = rgb(0,0,0,0.3)
colnoass = rgb(0,0,0,0.6)
par(mar = c(5,5,1,1)+0.1)
hist(log(Pnog[com10>0]),col=colass,main ='', ylab='Density', freq=FALSE, xlab='Log of probability',ylim=c(0, 0.45), xlim=c(-8,0), breaks=20, cex.lab=2, lty=1, lwd=4, cex.axis = 1.5)
hist(log(Pnog[com10==0]),col=colnoass,freq=FALSE, add=TRUE, lty=2, lwd=4)
legend(x='top', legend=c('Observed associations', 'Unobserved associations'), lwd=4, col=c(colass, colnoass),pt.cex=1.5, cex=1.5, lty=c(1,2))
dev.off()

## Degree Distribution with and without G
ZpostG = 1*(Pg>rocG$threshold)
ZpostNoG = 1*(Pnog>rocNoG$threshold)
pdf(paste0(subDir,'degree_dist.pdf'))
par(mfrow = c(2,2), mar = c(4.1, 4, 0.2, 0.2) +0.1)
plot_degree(com10, ZpostG, type='hosts',host.col = 'red')
plot_degree(com10, ZpostG, type='parasites', parasite.col = 'blue')
plot_degree(com10, ZpostNoG, type='hosts',host.col = 'red')
plot_degree(com10, ZpostNoG, type='parasites', parasite.col = 'blue')
dev.off()

    
## Printing top M interactions
## top m
m = 4000
## Full with G
ord.pg =order(Pg, decreasing=TRUE)
ord.pg[which(Pg[ord.pg]==1 & com10[ord.pg]==1)]
## Full without G
ord.p = order(Pnog, decreasing=TRUE)

topm = data.frame(t(sapply(1:m, function(r) c(actual = sum(com10[ord.pg[1:r]]),
    withG= sum(com10[ord.pg[1:r]]*ZpostG[ord.pg[1:r]]),
    withOutG= sum(com10[ord.p[1:r]]*ZpostNoG[ord.p[1:r]])))))

pdf(paste0(subDir, 'TopM_.pdf'))
par(mar = c(5, 5, 1, 0.2)+ 0.1)
plot(x=1:m,y = topm[,'withG'], xlab='Number of validated interactions', ylab = 'Number of recovered interactions',
     col='black',lty=1, type='l', lwd=3, cex.lab=2, cex.axis = 1.5, ylim = c(0,4000))
lines(x=1:m,y = topm[,'withOutG'], lty=5, type='l', lwd=3, col='black')
lines(x=1:m, y=1:m, lty=3, lwd=2, col='black')
gnames = c(paste('LS-net: ', MODEL, c('model with uncertainty', 'model')), 'x=y')
legend('topleft',legend = gnames,lty=c(1,5,3),lwd=c(3,3,2), cex=1.4)
dev.off()


## Saving workspace
if(SAVE_PARAM)
    save.image(file = paste0(subDir, SAVE_FILE))

##################################################
##################################################
