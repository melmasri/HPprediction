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
if(length(grep('com', ls()))==0)
    stop("no object named 'com' in the data file.")

if(length(grep('phy_dist', ls()))==0)
    stop("no object named 'phy_dist' in the data file.")

if(!is.matrix(com) | is.matrix(phy_dist))
    stop("either 'com' or 'phy_dist' are not a matrix type.")

if(!isSymmetric(phy_dist))
    stop("matrix 'phy_dist' is not symmetric.")

if(nrow(com)!= nrow(phy_dist))
    stop("matrix 'phy_dist' doesn't have the same number of rows as 'com'.")

## slice= ceiling(20000/ncol(com))
slice = SLICE
param = ICM_est(Z =1*(com>0),slice=slice ,dist= dist,eta=1,
    hyper = hyper, AdaptiveMC=TRUE,ICM.HORIZ= ICM.HORIZ)

if(SAVE_PARAM)
    save.image(file = 'param.RData')
##################################################
##################################################
