#########################################
## Script to run on the whole dataset
#rm(list= ls())
## Global Variable
SAVE_PARAM = TRUE
## DATAFILENAME = 'comGMPD.RData'
## DATAFILENAME = 'comEID-PS.RData'
## DATAFILENAME = 'comSim600x200.RData'
print(DATAFILENAME)

#source('library.R')
#source('gen.R')
load(DATAFILENAME)

#######################
## Parameters for the independent GGP
## set the correct prior.
print('Setting the prior.')
source('library')
## No uncertain
#A = generate_interactions(r=40, c=200, eta=1.5, aj=0.5, ai=0.7)
#plot_Z(A$Z)

slice = 100
param_phy = gibbs_one(Z=1*(com>0),slice=slice,dist =phy_dist, eta=1,uncertain=FALSE, hyper=hyper, AdaptiveMC= TRUE, updateHyper = FALSE)

if(SAVE_PARAM)
    save.image(file = 'param.RData')
##################################################
##################################################
