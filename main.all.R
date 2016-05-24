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
    hyper = list(parasite= c(1/3, 1), host =c(2, 1), eta = c(1.13, 1))

if(dataset =='eid')
    hyper = list(parasite= c(0.5, 1), host =c(0.1, 2), eta = c(0.96, 1))

#source('gen.R')
## Z = 1*(com>0)
## dist = Z%*%t(Z)
## dist  = 1-1/dist
## dist[is.infinite(dist)]<-0
## dist = dist + 1e-4
## diag(dist)<-0

## dist=phy_dist
## Summary statistics

yy = rgamma(200, shape=2, rate=1)
ww = rgamma(500, shape = 1/3, rate=1)
zz = 1*(runif(200*500)<= 1-exp(-outer(yy,ww)))
d = matrix(0, nrow=length(yy), ncol=length(yy))
d[upper.tri(d)]<-runif(length(yy)*(length(yy)-1)/2)
d = d + t(d)

aux = which(colSums(Z)==0)
zz<-zz[,-aux]
ww<-ww[-aux]
zz = 1-exp(-outer(yy,ww)*((d^2)%*%zz))
zz = 1*(runif(length(ww)*length(yy))<=zz)
plot_Z(zz)
com=zz
aux = which(colSums(Z)==0)
zz<-zz[,-aux]
aux = which(rowSums(Z)==0)
#zz<-zz[-aux,]
dim(zz)
com=zz
source('library.R')
source('gen.R')

param_phy = gibbs_one(Z=1*(com>0),slice=10, uncertain=FALSE, yMH=FALSE, wMH =FALSE, wEta = FALSE, yEta=FALSE)#, hyper= hyper)

if(SAVE_PARAM)
    save.image(file = 'param.RData')
##################################################
##################################################
