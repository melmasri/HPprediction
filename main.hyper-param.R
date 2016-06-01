## Global Variable
args <- commandArgs(trailingOnly = TRUE)
print(args)

dataset= args[1]

SAVE_PARAM = TRUE

if(dataset=='gmp')
    DATAFILENAME = 'comGMPD.single.RData'
if(dataset == 'eid')
    DATAFILENAME = 'comEID-PS..single.RData'


sink(paste0('Dataset: ', dataset))
print(DATAFILENAME)
source('library.R')
source('gen.R')
load(DATAFILENAME)

slice= min(5,ceiling(10000/ncol(com)))

if(dataset =='gmp')
    hyper = list(parasite =c(5, 1), host = c(1,1), eta = c(0.008)) 

if(dataset =='eid')
    hyper = list(parasite= c(5, 1), host =c(1, 2), eta = c(0.01))

param_phy1 = gibbs_one(Z=1*(com>0),slice=slice,dist = phy_dist, eta=1,
    uncertain=FALSE, yMH=FALSE, wMH =FALSE, wEta = FALSE, yEta=FALSE, hyper=hyper,
    updateHyper=TRUE, AdaptiveMC=TRUE)

param_phy1[['w']]<-NULL
param_phy1[['y']]<-NULL
param_phy1[['eta']]<-NULL
param_phy1[['g']]<-NULL

param_phy = list(w = w0, y = y0, burn_in = burn_in - max(-throw.out), throw.out = max(-throw.out),eta = eta, g=g, hh=hh, sd = list(w=w_sd, y = y_sd, eta= eta_sd))

if(dataset =='gmp')
    hyper = list(parasite =c(1, 1), host = c(2,1), eta = c(0.008)) # comGMPD.single.RData

if(dataset =='eid')
    hyper = list(parasite= c(1, 1), host =c(2, 2), eta = c(0.01))

param_phy2 = gibbs_one(Z=1*(com>0),slice=slice,dist = phy_dist, eta=1,
    uncertain=FALSE, yMH=FALSE, wMH =FALSE, wEta = FALSE, yEta=FALSE, hyper=hyper,
    updateHyper=TRUE, AdaptiveMC=TRUE)

param_phy2[['w']]<-NULL
param_phy2[['y']]<-NULL
param_phy2[['eta']]<-NULL
param_phy2[['g']]<-NULL

if(dataset =='gmp')
    hyper = list(parasite =c(0.1, 1), host = c(0.1,1), eta = c(0.008)) 

if(dataset =='eid')
    hyper = list(parasite= c(0.1, 1), host =c(0.1, 2), eta = c(0.01))

param_phy3 = gibbs_one(Z=1*(com>0),slice=slice,dist = phy_dist, eta=1,
    uncertain=FALSE, yMH=FALSE, wMH =FALSE, wEta = FALSE, yEta=FALSE, hyper=hyper,
    updateHyper=TRUE, AdaptiveMC=TRUE)

param_phy3[['w']]<-NULL
param_phy3[['y']]<-NULL
param_phy3[['eta']]<-NULL
param_phy3[['g']]<-NULL

if(SAVE_PARAM)
    save.image(file = paste0('param-',dataset,'.RData'))

##################################################
##################################################
sink()
system("grep '^[^>+;]' report.txt > report_clean.txt")
setwd('../')
subj = '"End of sim "'
body = 'End HyperChain'
email = paste("echo '" ,body,"' | mailx -s ", subj, " elmasri.m@gmail.com")
system(email)
q('no')
