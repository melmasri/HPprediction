#########################################
## Script to run on the whole dataset
rm(list= ls())

## Loading required packages
library(ape)
library(geiger)
library(fulltext)

## loading mammal supertree included in Fritz et al. (2009)
source('example-GMPD/download_tree.R')       # see variable 'tree'
pruned.tree <- drop.tip(tree,
                        sample(tree$tip.label)[1:(0.8*length(tree$tip.label))]) # trimming 80% of tree tips for speed

## loading GMPD
source('example-GMPD/load_GMPD.R')           # see matrix 'com'

## sourcing MCMC script
source('network_MCMC.R')

## running the model of interest
obj = network_est(Z = com, slices=1000, tree=pruned.tree, model.type='full') # full model
names(obj)
names(obj$param)

## load useful network analysis functions
source('network_analysis.R')

## plot input matrix and tree
plot(obj$tree, show.tip.label=FALSE)
plot_Z(obj$Z)

## Probability matrix
## Extracting mean posteriors
Y = if(is.matrix(obj$param$y)) rowMeans(obj$param$y) else  mean(obj$param$y)
W = if(is.matrix(obj$param$w)) rowMeans(obj$param$w) else  mean(obj$param$w)
Eta = if(is.null(obj$param$eta)) 0 else mean(obj$param$eta)

## Affinity only model
P =  1-exp(-outer(Y, W))

### Full or distance model
## Creating distance
distance = 1/cophenetic(rescale(obj$tree, 'EB', Eta))
diag(distance)<-0
distance = distance %*% obj$Z
distance[distance==0]<-1

## Probability matrix
## Full model
P = 1-  exp(-outer(Y, W)*distance)

## Distance only model 
P = 1-  exp(-distance)

## ROC curves, AUC and posterior Z
roc = rocCurves(obj$Z, obj$Z, P = P, all = TRUE) # ROC
plot_Z(P>roc$threshold + 0)                      # posterior Z

##################################################
##################################################
