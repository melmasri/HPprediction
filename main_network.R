#########################################
## Script to run on the whole dataset
rm(list= ls())
library(ape)
library(geiger)


## loading tree
source('example/download_tree.R')

## loading data
source('example/load_GMPD.R')

## sourcing MCMC script
source('networkMCMC.R')


obj = network_est(Z = com, slices=5, tree=tree, model.type='dist', burn.in=0)

## Probability matrix
## Extracting mean posteriors
Y = ifelse(is.matrix(obj$param$y), rowMeans(obj$param$y), mean(obj$param$y))
W = ifelse(is.matrix(obj$param$w), rowMeans(obj$param$w), mean(obj$param$w))
Eta = mean(obj$param$eta)

## Affinity only model
P =  1-exp(-outer(Y, W))

### Full model (both) or Distance only model
## Creating distance
distance = 1/cophenetic(rescale(obj$tree, 'EB', Eta))
diag(distance)<-0
distance = distance %*% obj$Z
distance[distance==0]<-1

## Probability matrix
## Full model (both)
P = 1-  exp(-outer(Y, W)*distance)

## Distance only model 
P = 1-  exp(-distance)


## Analysis
source('library.R')

## rocCurves 
roc = rocCurves(obj$Z, obj$Z, P = P, all = TRUE)


##################################################
##################################################
