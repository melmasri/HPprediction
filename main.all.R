#########################################
## Script to run on the whole dataset
#rm(list= ls())
## Global Variable
SAVE_PARAM = TRUE
## DATAFILENAME = 'comGMPD.RData'
## DATAFILENAME = 'comEID-PS.RData'
print(DATAFILENAME)

## source('library.R')
## source('gen.R')
load(DATAFILENAME)
#######################
## print tests
library(ape)
library(geiger)

if(length(grep('com', ls()))==0)
    stop("no object named 'com' in the data file.")

if(length(grep('tree', ls()))==0)
    stop("no object named 'tree' in the data file.")

tree <- read.tree('../mammals.tre')
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(com)])

# Testing all names in hosts in com exist in tree
if(! all(sapply(rownames(com), function(r) r %in% tree$tip.label))) {
    print('Warning! Not all hosts in com exist in tree. Hosts not found in tree will be removed.')
    com <- com[rownames(com)%in%tree$tip.label,]
}

dd = cophenetic(rescale(tree, 'EB', 0))
host.order <- sapply(rownames(dd), function(r) which(r==rownames(com)))
com = com[host.order,]


param = ICM_est(Z =1*(com>0),slice=SLICE, tree=tree, eta=0, beta=0.5,
    eta_sd=0.005, a_y = 0.15, a_w= 0.15)
    
if(SAVE_PARAM)
    save.image(file = 'param.RData')
##################################################
##################################################
