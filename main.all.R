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

if(dataset =='gmp')
    hyper = list(parasite =c(29.8, 1), host = c(0.24,1), eta = c(0.008)) # comGMPD.single.RData

if(dataset =='eid')
    hyper = list(parasite= c(0.5, 1), host =c(0.1, 2), eta = c(0.01))

#source('gen.R')
## Z = 1*(com>0)
## dist = Z%*%t(Z)
## dist  = 1-1/dist
## dist[is.infinite(dist)]<-0
## dist = dist + 1e-4
## diag(dist)<-0

## dist=phy_dist
## Summary statistics


## source('library.R')
## source('gen.R')

## y = rgamma(100, shape=2, rate=3)
## w = rgamma(500, shape=0.5, rate=0.5)
## z = 1*(runif(100*500)<=1-exp(-outer(y,w)))
## com=z
## SIMPLERHO=TRUE

param_phy = gibbs_one(Z=1*(com>0),slice=8,dist = phy_dist, eta=1,uncertain=FALSE, yMH=FALSE, wMH =!SIMPLERHO, wEta = !SIMPLERHO, yEta=FALSE, hyper=hyper, updateHyper=FALSE)

if(SAVE_PARAM)
    save.image(file = 'param.RData')
##################################################
##################################################
