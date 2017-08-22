#########################################
## Script to run on the whole dataset
rm(list= ls())
library(ape)
library(geiger)


## loading tree
source('example/download_tree.R')

## loading data
source()

## sourcing MCMC script
source('networkMCMC.R')

param = network_est(Z = com, slices=5, tree=tree, model.type='both')




##################################################
##################################################
