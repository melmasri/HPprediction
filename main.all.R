#########################################
## Script to run on the whole dataset
#rm(list= ls())
## Global Variable
SAVE_PARAM = TRUE
## DATAFILENAME = 'comGMPD.RData'
## DATAFILENAME = 'comEID-PS.RData'
print(DATAFILENAME)

#source('library.R')
#source('gen.R')
load(DATAFILENAME)

#######################
## Parameters for the independent GGP
## set the correct prior.
print('Setting the prior.')
## No uncertain

## if(dataset =='gmp')
##     hyper = list(parasite =c(29.8, 1), host = c(0.24,1), eta = c(0.008)) # comGMPD.single.RData

## if(dataset =='eid')
##     hyper = list(parasite= c(0.5, 1), host =c(0.1, 2), eta = c(0.01))

slice= ceiling(20000/ncol(com))
param_phy = gibbs_one(Z=1*(com>0),slice=slice,dist = phy_dist, eta=1,uncertain=FALSE, hyper=hyper, AdaptiveMC= TRUE, updateHyper = FALSE)

if(SAVE_PARAM)
    save.image(file = 'param.RData')
##################################################
##################################################
