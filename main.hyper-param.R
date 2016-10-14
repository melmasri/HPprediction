## Global Variable
## TO run
## Rscript main.hyper-param.R dataset
#dataset: gmp or eid

args <- commandArgs(trailingOnly = TRUE)
print(args)

dataset= args[1]

SAVE_PARAM = TRUE

if(dataset=='gmp')
    DATAFILENAME = 'comGMPD.single.RData'
if(dataset == 'eid')
    DATAFILENAME = 'comEID-PS.single.RData'

subDir = paste0('hyper-',dataset)
## Creating a subdirectory
dir.create(file.path(subDir))
## Setting the working directory
setwd(file.path(subDir))
## starting the process

sink(paste0('Dataset: ', dataset))
print(DATAFILENAME)
source('../library.R')
source('../gen.R')
load(DATAFILENAME)

slice= min(8,ceiling(4000/ncol(com)))

if(dataset =='gmp')
    hyper = list(parasite =c(5, 1), host = c(1,1), eta = c(0.008)) 

if(dataset =='eid')
    hyper = list(parasite= c(5, 1), host =c(1, 2), eta = c(0.01))

print('hyper settings 1:')
print(hyper)

param_phy1 = gibbs_one(Z=1*(com>0),slice=slice,dist = phy_dist, eta=1,
    uncertain=FALSE, hyper=hyper,updateHyper=TRUE, AdaptiveMC=TRUE)

param_phy1[['w']]<-NULL
param_phy1[['y']]<-NULL
param_phy1[['eta']]<-NULL
param_phy1[['g']]<-NULL

if(dataset =='gmp')
    hyper = list(parasite =c(1, 1), host = c(0.1,1), eta = c(0.008)) # comGMPD.single.RData

if(dataset =='eid')
    hyper = list(parasite= c(1, 1), host =c(2, 2), eta = c(0.01))

print('hyper settings 2:')
print(hyper)

param_phy2 = gibbs_one(Z=1*(com>0),slice=slice,dist = phy_dist, eta=1,
    uncertain=FALSE, hyper=hyper, updateHyper=TRUE, AdaptiveMC=TRUE)

param_phy2[['w']]<-NULL
param_phy2[['y']]<-NULL
param_phy2[['eta']]<-NULL
param_phy2[['g']]<-NULL

if(dataset =='gmp')
    hyper = list(parasite =c(0.1, 1), host = c(0.3,1), eta = c(0.008)) 

if(dataset =='eid')
    hyper = list(parasite= c(0.1, 1), host =c(0.1, 2), eta = c(0.01))

print('hyper settings 2:')
print(hyper)

param_phy3 = gibbs_one(Z=1*(com>0),slice=slice,dist = phy_dist, eta=1,
    uncertain=FALSE,hyper=hyper, updateHyper=TRUE, AdaptiveMC=TRUE)

param_phy3[['w']]<-NULL
param_phy3[['y']]<-NULL
param_phy3[['eta']]<-NULL
param_phy3[['g']]<-NULL

if(SAVE_PARAM)
    save.image(file = paste0('param-',dataset,'.RData'))



par(mfrow=c(2,1))

pdf('hyper-param-hosts.pdf')
plot(param_phy$hh[1,], type='l', main = '', col='blue', xlab = 'Iteration')
lines(param_phy2$hh[1,], col='orange')
lines(param_phy3$hh[1,],col ='yellow')
dev.off()

pdf('hyper-param-hosts.pdf')
plot(param_phy$hh[3,], type='l', main = '', col='blue', xlab = 'Iteration')
lines(param_phy2$hh[3,], col='orange')
lines(param_phy3$hh[3,],col ='yellow')
dev.off()

burn =floor(param_phy$burn_in/2):param_phy$burn_in

a = rbind(rowMeans(param_phy$hh[,burn]), rowMeans(param_phy2$hh[,burn]),rowMeans(param_phy3$hh[,burn]))

colnames(a)<-c('host-alpha', 'host-tau', 'parasite-alpha', 'parasite-tau')
print('Hyper parameters:')
print(a)

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
