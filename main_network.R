#########################################
## Script to run on the whole dataset
#########################################

## General variables
## please specify the following parameters
## SAVE_PARAM = TRUE                    # should workspace be saved
## SAVE_FILE = 'param.RData'            # name of output R workspace file
## MODEL = 'full'                       # full, distance or affinity
## SLICE = 100                          # no of iterations
## subDir = ''                          # directory to print the results 
## COUNT = TRUE                         # TRUE = count data, FALSE = year of first pub.

## loading mammal supertree included in Fritz et al. (2009)
source('example-GMPD/download_tree.R')       # see variable 'tree'

## loading GMPD
source('example-GMPD/load_GMPD.R')           # see matrix 'com'

## sourcing MCMC script
source('network_MCMC.R')

## running the model of interest
obj = network_est(Z = com, slices=SLICE, tree=tree, model.type=MODEL, backup=TRUE) # full model
names(obj)
names(obj$param)

## load useful network analysis functions
source('network_analysis.R')

## Probability matrix
## Extracting mean posteriors
Y = if(is.matrix(obj$param$y)) rowMeans(obj$param$y) else  mean(obj$param$y)
W = if(is.matrix(obj$param$w)) rowMeans(obj$param$w) else  mean(obj$param$w)
Eta = if(is.null(obj$param$eta)) 0 else mean(obj$param$eta)
YW = outer(Y, W)
## setting affinity to 1 in distance model
if(grepl('dist', MODEL)) YW = 1 

### Full or distance model
## Creating distance
if(grepl('(full|dist)', MODEL)){
    distance = 1/cophenetic(rescale(obj$tree, 'EB', Eta))
    diag(distance)<-0
    distance = distance %*% obj$Z
    distance[distance==0] <- if(grepl('dist', MODEL)) Inf else 1
} else distance = 1

## models
## P = 1-exp(-outer(Y, W))                 # affinity model
## P = 1-exp(-distance)                    # distance model
## P = 1-exp(-YW*distance)                 # full model

## All in one probability matrix
P = 1-  exp(-YW*distance)

## ROC curves and AUC
roc = rocCurves(obj$Z, obj$Z, P = P, all = TRUE, plot = FALSE) # ROC

## some numerical analysis
TB = ana.table(obj$Z, obj$Z, P, roc, TRUE)
## Printing and writing out average MCMC 
print(sprintf('Model: %s, AUC: %f and percent 1 recovered out of all: %f',
              MODEL,mean(TB$auc), mean(TB$pred.tot.ones)))

write.csv(TB, file = paste0(subDir, 'AUC-PRED.csv'))


## left ordering and outputting interaction matrix
indices = lof(obj$Z, indices = TRUE)
## printing input interaction matrix
pdf(paste0(subDir, 'Z_input.pdf'))
plot_Z(obj$Z[, indices])
dev.off()
## printing input tree
pdf(paste0(subDir, 'tree_input.pdf'))
plot(obj$tree, show.tip.label=FALSE)
dev.off()
## printing posterior interaction matrix
pdf(paste0(subDir, 'Z_', MODEL, '.pdf'))
plot_Z(1*(P[, indices]>roc$threshold))                     
dev.off()




## MCMC trace plots and ACF
require(coda)
pdf(paste0(subDir, 'param_mcmc_acf.pdf'))
r = which.max(rowSums(com))
c = which.max(colSums(com))
if(grepl('aff', MODEL)) par(mfrow=c(2,2))
if(grepl('dist', MODEL)) par(mfrow=c(1,2))
if(grepl('full', MODEL)) par(mfrow=c(3,2))
if(grepl('(aff|full)', MODEL)){
    plot(obj$param$y[r,], type='l', main='',xlab='Iteration', ylab='Host')
    plot(obj$param$w[c,], type='l', main='', xlab = 'Iteration', ylab='Parasite')
}
if(grepl('(dist|full)', MODEL))
    plot(obj$param$eta, type='l', main = '', xlab = 'Iteration', ylab='Tree scaling parameter')
if(grepl('(aff|full)', MODEL)){
    acf(obj$param$y[r,],lag.max=100, main='Host - most active')
    text(80, 0.8, paste0('Effective sample size: ',round(effectiveSize(obj$param$y[r,]))))
    round(effectiveSize(obj$param$y[r,]))
    acf(obj$param$w[c,],lag.max=100, main='Parasite - most active')
    text(80, 0.8, paste0('Effective sample size: ',round(effectiveSize(obj$param$w[c,]))))
    round(effectiveSize(obj$param$w[c,]))
}
if(grepl('(dist|full)', MODEL)){
    acf(obj$param$eta,lag.max=100, main='Scale parameter')
    text(80, 0.8, paste0('Effective sample size: ',round(effectiveSize(obj$param$eta))))
    round(effectiveSize(obj$param$eta))
}
dev.off()


## Saving workspace
if(SAVE_PARAM)
    save.image(file = paste0(subDir, SAVE_FILE))

##################################################
##################################################


