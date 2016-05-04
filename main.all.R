# ########################################################################
# Mohamad Elmasri (elmasri.m@gmail.com)
# This is script is supposed to replicate the results of:
#  Caron, F. (2012) Bayesian nonparametric models for bipartite graphics
# ########################################################################
# Project dates:	start January 19, 2015
# 					close ongoing
#########################################
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
if(dataset =='gmp'){
    a_w = 0.6
    b_w = 1
    a_y = 0.5
    b_y = 1
    a_e = 1.13
    b_e = 1
    eta_sd = 0.006
    w_sd = 0.15
    y_sd = 0.15
}
if(dataset =='eid'){
    a_w = 0.5
    b_w = 1
    a_y = 0.1
    b_y = 2
    a_e = 0.96
    b_e = 1
    eta_sd = 0.003
    w_sd = 0.15
    y_sd = 0.15
}
#source('gen.R')
## Z = 1*(com>0)
## dist = Z%*%t(Z)
## dist  = 1-1/dist
## dist[is.infinite(dist)]<-0
## dist = dist + 1e-4
## diag(dist)<-0

## dist=phy_dist
## Summary statistics
param_phy = gibbs_one(1*(com>0),slice=10,dist= phy_dist, eta=1, uncertain=FALSE, yMH=TRUE, wMH = TRUE)

if(SAVE_PARAM)
    save.image(file = 'param.RData')

##################################################
##################################################
