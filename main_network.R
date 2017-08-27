#########################################
## Script to run on the whole dataset
rm(list= ls())

## Loading required packages
library(ape)
library(geiger)
library(fulltext)


## loading mammal supertree included in Fritz et al. 2009 (DOI: 10.1111/j.1461-0248.2009.01307.x). 
source('example/download_tree.R')       # see variable tree

## loading GMPD
source('example/load_GMPD.R')           # see matrix com

## sourcing MCMC script
source('networkMCMC.R')


obj = network_est(Z = com, slices=100, tree=tree, model.type='affinit')

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


## Analysis
source('library.R')

## rocCurves 
roc = rocCurves(obj$Z, obj$Z, P = P, all = TRUE)


##################################################
##################################################
