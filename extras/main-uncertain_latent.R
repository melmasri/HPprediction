#########################################
## Script to run using cross-validation 
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

## Script for the nearest NN

MODEL = 'latent-uncertain'
COUNT = FALSE
SINGLE = TRUE
RAND  = TRUE
SAVE_PARAM = TRUE
SAVE_FILE = 'param.RData'
sTime = Sys.time()
## Formating the sub-directory name
subDir =  paste(toupper('cv'),toupper(MODEL), format(sTime, "%d-%m-%Hh%M"), sep='-')

## adding sub-extension
subDir = paste0('./', subDir, '/')

dir.create(file.path(subDir))
## report file
reportFile <- paste0(subDir, "report.txt")
## Setting the working directory

## starting the process
sink(reportFile)
print(date())

set.seed(23456)
## Loading required packages
library(parallel)
library(latentnet)
## loading mammal supertree included in Fritz et al. (2009)
source('example-GMPD/download_tree.R')  # see variable 'tree'

## loading GMPD
source('example-GMPD/load_GMPD.R')      # see matrix 'com'
if(!SINGLE){
    aux = which(colSums(1*(com>0))==1)
    com = com[, -aux]
    com = com[-which(rowSums(1*(com>0))==0), ]
}
dim(com)
## sourcing MCMC script
source('network_MCMC.R')

## preparing tree and com
cleaned = network_clean(com, tree, 'full')
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

folds = if(RAND)
            cross.validate.fold(com, n= 5, 1, 'random') else
cross.validate.fold(com, n= 5, 1) 
                                        
tot.gr = length(unique(folds[,'gr']))   # total number of CV groups

## A loop ran over all CV groups
res = mclapply(1:tot.gr ,function(x, folds, Z, tree, slice, model.type){
    ## Analysis for a single fold
    Z.train = Z
    Z.train[folds[which(folds[,'gr']==x),c('row', 'col')]]<-0
    ## running the model of interest
    X = network(Z.train)
    fit<-ergmm(X~euclidean(d=2)+ rsociality,
               control=ergmm.control(mle.maxit=10,burnin=0,threads=2),verbose=TRUE)
    pred <- predict(fit)

    parasites = which(network.vertex.names(X) %in% colnames(Z))
    hosts = which(network.vertex.names(X) %in% rownames(Z))
    P = matrix(0, nrow(Z), ncol(Z))
    P= t(pred[parasites, hosts])
    
    ## order the rows in Z.test as in Z.train
    list(param=P)
},folds=folds,Z = com, tree=tree, model.type=MODEL, slice = SLICE,
    mc.preschedule = TRUE, mc.cores = min(tot.gr, NO.CORES))

## Constructing the P probability matrix from CV results
aux = rowMeans(sapply(res, function(r) r$param))
P = matrix(aux, nrow(com), ncol(com))

## ROC curves and posterior interaction matrices
roc = rocCurves(com10, com,P, all=TRUE, plot=FALSE, bins=400)
pdf(paste0(subDir, 'latent_Z-2010.pdf'))
plot_Z(1*(P> roc$threshold))
dev.off()

## Printing AUC and prediction results
tb = ana.table(com10, com, P,   roc = roc)
print(tb)

if(SAVE_PARAM)
    save.image(file = paste0(subDir, SAVE_FILE))


##################################################
##################################################
##-----------------------------------------------------------
eTime = Sys.time()
print(sprintf('End time %s.',format(eTime, "%Y-%m-%d %H:%M:%S")))
## Processing time
print('Total processing time')
print(eTime - sTime)

######################
## Closing sink and reverting work directory.
sink()
system(paste("grep '^[^>+;]'", reportFile, ">", paste0(subDir, "report_clean.txt") ))


q('no')
