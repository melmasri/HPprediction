# ########################################################################
# Mohamad Elmasri (elmasri.m@gmail.com)
# This is script is supposed to replicate the results of:
#  Caron, F. (2012) Bayesian nonparametric models for bipartite graphics
# ########################################################################
# Project dates:	start January 19, 2015
# 					close ongoing
#########################################
## rm(list= ls())
## Global Variable
SAVE_PARAM = TRUE
PDFPLOT = TRUE                          # If user want to save to
DATAFILENAME = 'comEID-PS.RData'
print(DATAFILENAME)
source('library.R')
source('gen.R')
load(DATAFILENAME)
com_pa = 1*(com>0)
com_paCross = 1*(comCross>0)
#pairs = cross.validate.fold(com)

#######################
## Parameters for the independent GGP
## set the correct prior.
print('Setting the prior.')
if(dataset =='gmp'){
    a_w = 0.6
    b_w = 1
    a_y = 0.5
    b_y = 1
    a_e = 1.13
    b_e = 1
}else{
    a_w = 0.5
    b_w = 1
    a_y = 0.1
    b_y = 2
    a_e = 0.96
    b_e = 1
}
print(sprintf('a_w %0.3f, b_w %0.3f, a_y %0.3f, b_y %0.3f, a_e %0.3f, b_e %0.3f',a_w, b_w, a_y, b_y, a_e, b_e))

                                        #sigma = 0
#tau =1
#a = 8
#process ='gamma'	 		#Allowed c('gamma','gen.gamma','inv.gamma','stable')
								# gamma 	(sigma=0)
                                # inv.gamma (sigma=0.5)
								# gen.gamma (sigma>0)
                                # stable 	(tau=0)

## Generative process as in Section 2.4

## nodes_w 	= ncol(com_pa)		    # Number of books
## nodes_y		= nrow(com_pa)			# Number of readers
#a_y =0.2;b_y =tau ;
## y <-rep(2,nodes_y);
##y<-rgamma(nodes_y,a_y,b_y)		 # interest parameter
## y <- rep(2, nodes_y)

## Estimation
##-----------------------------------------------------------
## Summary statistics
cnames = colnames(comCross)
rnames = rownames(comCross)
colnames(comCross)<-1:ncol(comCross)
rownames(comCross)<-1:nrow(comCross)
colnames(phy_dist)<-1:ncol(phy_dist)
rownames(phy_dist)<-1:nrow(phy_dist)

param_phy1 = gibbs_one(com_pa,slice=17 ,dist= phy_dist, eta=1, uncertain=FALSE, yMH=FALSE, wMH = TRUE)
param_phy = param_phy1

## Analysis
r = which.max(rowSums(com_paCross));r
c = which.max(colSums(com_paCross));c

if(PDFPLOT) pdf(paste0(DATAFILENAME, 'plot.pdf'))

burn = -c(1:2)
par(mfrow=c(3,1))
plot(param_phy$y[r,burn], type='l', main = 'Most popular host')
plot(param_phy$w[c,burn], type='l', main = 'Most popular parasite')
if(!is.null(param_phy$eta))
    plot(param_phy$eta[burn], type='l', main = 'Beta posterior')

par(mfrow=c(2,1))
if(!is.null(param_phy$L)){
    plot(param_phy$L[1,], type='l', main = 'L1 posterior')
    plot(param_phy$L[2,], type='l', main = 'L0posterior')
}

par(mfrow=c(1,2))
if(!is.null(param_phy$Z) & !is.null(param_phy$L)){
    hist(param_phy$L['l1',] - g.old, xlim = c(0,max(param_phy$L['l1',])))
    hist(param_phy$L['l0',], add=T, col='red')
    plot(density(param_phy$L['l1',]))
}

##P = get.P.mode(param_phy,com_paCross, phy_dist, lambda_phy)
print('Regular trans aux')
param_phy$w = param_phy$w^(1/param_phy$eta)
## aux= getMode(param_phy)
aux = getMean(param_phy)
print("aux$eta");print(aux$eta)
P1 = 1-  exp(-outer(aux$y, aux$w^aux$eta)*((phy_dist^aux$eta)%*% com_paCross))

roc = rocCurves(Z =com_pa, Z_cross = com_paCross, P=P1)
plot_Z(1*(P1>roc$threshold), xlab= paste('Threshold: ', round(roc$threshold,3)))
if(PDFPLOT) dev.off()

zz = 1*(P1>roc$threshold)
zz1 = sum(zz[com_pa==1 & com_paCross==0])/sum(abs(com_paCross - com_pa)[com_pa==1])

samples.table = data.frame(pa =c(sum(zz - com_pa ==1),  sum(zz - com_pa == -1)),
   paCA =  c(sum(zz - com_paCross ==1), sum(zz - com_paCross ==-1) ))
rownames(samples.table)<-c('new.in', 'old.out ')
print(samples.table)

print(sprintf('Total interactoins %i, hold out is %i samples of which %0.3f predicted',sum(com_pa),sum(abs(com_paCross - com_pa)[com_pa==1]),zz1))
print(sprintf('AUC is %0.2f and probability threshold is %0.3f', roc$auc, roc$threshold))

x = topPairs(P=P1, Z= com_paCross, topX=40)
zz = melt(com_pa)
x= cbind(x,ture = zz[as.numeric(rownames(x)),'value'])
print('Top 40 candidates and their true value')
print(x)
if(SAVE_PARAM)
    save.image(file = 'param.RData')

##################################################
##################################################
