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
if(exists("PATH.TO.FILE") && !is.null(PATH.TO.FILE)){
    if(grepl('.rds', PATH.TO.FILE, ignore.case = TRUE))
        com <- readRDS(PATH.TO.FILE) else 
    load(PATH.TO.FILE)
}else{
    source('example-GMPD/load_GMPD.R')           # see matrix 'com'    
}


## sourcing MCMC script
source('network_MCMC.R')

## running the model of interest
obj = network_est(Z = com, slices=SLICE, tree=tree, model.type=MODEL,
                  a_y=ALPHA.ROWS, a_w= ALPHA.COLS, burn.in=0) # full model
names(obj)
names(obj$param)

## load useful network analysis functions
source('network_analysis.R')

## Probability matrix
## Extracting mean posteriors of P
P = sample_parameter(obj$param, MODEL, obj$Z, obj$tree)
    

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

## printing output tree
if(grepl('(full|dist)', MODEL)){
    Eta = mean(obj$param$eta)
    pdf(paste0(subDir, 'tree_', MODEL,'.pdf'))
    plot(rescale(obj$tree, 'EB', Eta), show.tip.label=FALSE)
    dev.off()
}

## printing posterior interaction matrix
pdf(paste0(subDir, 'Z_', MODEL, 'new.pdf'))
plot_Z(1*(P[, indices]>roc$threshold))                     
dev.off()

## MCMC trace plots and ACF
require(coda)
pdf(paste0(subDir, 'param_mcmc_acf.pdf'))
r = which.max(rowSums(com))
c = which.max(colSums(com))
par(mar = c(5,5,1,1)+0.1)
if(grepl('aff', MODEL)) par(mfcol=c(2,2))
if(grepl('dist', MODEL)) par(mfcol=c(1,2))
if(grepl('full', MODEL)) par(mfcol=c(3,2))
if(grepl('(aff|full)', MODEL)){
    plot(obj$param$y[r,], type='l', main='',xlab='Iteration', ylab='Host',  cex.lab=2, cex.axis=1.5)
    plot(obj$param$w[c,], type='l', main='', xlab = 'Iteration', ylab='Parasite',  cex.lab=2, cex.axis=1.5)
}
if(grepl('(dist|full)', MODEL))
    plot(obj$param$eta, type='l', main = '', xlab = 'Iteration', ylab='Scale', cex.axis=1.5, cex.lab=2)
if(grepl('(aff|full)', MODEL)){
    acf(obj$param$y[r,],lag.max=100, main='', cex.axis=1.5, cex.lab=1.5)
    text(50, 0.9, paste0('Effective sample size: ',round(effectiveSize(obj$param$y[r,]))), cex = 1.5)
    round(effectiveSize(obj$param$y[r,]))
    acf(obj$param$w[c,],lag.max=100, main='', cex.axis=1.5, cex.lab=1.5)
    text(50, 0.9, paste0('Effective sample size: ',round(effectiveSize(obj$param$w[c,]))), cex = 1.5)
    round(effectiveSize(obj$param$w[c,]))
}
if(grepl('(dist|full)', MODEL)){
    acf(obj$param$eta,lag.max=100, main='', cex.axis=1.5, cex.lab=1.5)
    text(50, 0.9, paste0('Effective sample size: ',round(effectiveSize(obj$param$eta))), cex = 1.5)
    round(effectiveSize(obj$param$eta))
}
dev.off()


## Saving workspace
if(SAVE_PARAM)
    save.image(file = paste0(subDir, SAVE_FILE))

##################################################
##################################################


