#########################################
## Script to run on the whole dataset
rm(list= ls())
library(ape)
library(geiger)


## loading tree
source('example/download_tree.R')

## source('library.R')
## source('gen.R')


tree <- read.tree('../Data/mammals.tre')
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(com)])



param = network_est(Z = com, slices=5, tree=tree, model.type='both')

    
if(SAVE_PARAM)
    save.image(file = 'param.RData')
##################################################
##################################################
aux = getMean(param$param)
pdist = 1/cophenetic(rescale(param$input$tree, 'EB', aux$eta))
diag(pdist)<-0
db = pdist%*%param$input$Z
P = 1- exp(-outer(aux$y, aux$w)*db)

ana.plot(param$param)
roc = rocCurves(param$input$Z,param$input$Z,plot=TRUE,P,bins=400, all=TRUE)

    dim(outer(aux$y, aux$w))
dim(cophenetic(rescale(param$input$tree, 'EB', aux$eta)))
plot_Z(1*(P>roc$thres))
