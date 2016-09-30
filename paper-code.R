
## Plotting Hist of counts
##GMP
library(xtable)
source('library.R')

com0 = com[com>0]
r = range(com)
h = hist(com0, xlab='count',ylab='', main='',freq=TRUE, breaks =r[2])

pdf('GMP_com_hist_log.pdf')
plot(h$breaks[-r[2]],h$count, log="y", type='h', lwd=10, lend=2, col='gray', xlab='count', ylab ='log frequency')
dev.off()

##EID
## All counte above 704 have one association.
com0 = com[com>0]
x = table(com0)
names(x)[x>1]
com0 = com[com<=705]
r = range(com0)
h = hist(com0, xlab='count',ylab='', main='',freq=TRUE, breaks =r[2])

pdf('EID_com_hist_log.pdf')
plot(h$breaks[-r[2]],h$count, log="y", type='h', lwd=10, lend=2, col='gray', xlab='count', ylab ='log frequency')
dev.off()


## Plottig Z and degree
## GMP
load('comGMPD.RData')
source('library.R')

com_pa<-com
com_pa[com_pa>1]<-1
# Z
com_lof <-lof(com)
pdf('GMP_Z.pdf')
plot_Z(com_lof>0 + 0, xlab = 'parasites', ylab = 'hosts')
dev.off()

pdf('GMP_degree.pdf')
plot_degree(com_lof>0 + 0)
dev.off()

## EID
load('comEID.RData')
source('library.R')

com_pa<-com
com_pa[com_pa>1]<-1
# Z
com_lof <-lof(com)
pdf('EID_Z.pdf')
plot_Z(com_lof>0 + 0, xlab = 'parasites', ylab = 'hosts')
dev.off()

pdf('EID_degree.pdf')
plot_degree(com_lof>0 + 0)
dev.off()


##################################################
##################################################
##################################################
### Section: paramter estimation
rm(list=ls())
## load('gmp-11-03-08h43/param.RData')
## load('gmp-30-04-21h50/param.RData')
## load('gmp-30-04-23h34/param.RData')
load('gmp-01-05-00h46/param.RData')
## load('eid-11-03-08h43/param.RData')
## load('eid-30-04-21h49/param.RData')
load('eid-01-05-00h46/param.RData')
load('/home/mo/Github/HP-prediction/gmp-24-05-14h48/param.RData')
library(xtable)

com_pa = 1*(com>0)
r = which.max(rowSums(com_pa));r
c = which.max(colSums(com_pa));c
burn = param_phy$burn_in - 20000:1
dim(param_phy$y)
burn = burn[burn>0]
range(burn)

pdf(paste0(dataset, '_Z.pdf'))
plot_Z(1*(com>0) , xlab = 'parasites', ylab = 'hosts')
dev.off()

pdf(paste0(dataset, '_degree.pdf'))
plot_degree(com>0 + 0)
dev.off()


pdf(paste0(dataset, '_param_mcmc.pdf'))

par(mfrow=c(3,1))
plot(param_phy$y[r,burn], type='l', main='', col='blue', xlab='Iteration', ylab='Host parameter')
plot(param_phy$w[c,burn], type='l', main = '', col='blue', xlab = 'Iteration', ylab='Parasite parameter')
plot(param_phy$eta[burn], type='l', main = '', col='blue', xlab = 'Iteration', ylab='Scaling parameter')

## plot(param_phy$y[r,burn]*param_phy$y[c, burn], type='l', main='', col='blue', xlab='Iteration', ylab='Host parameter')

par(mfrow=c(2,1))
plot(param_phy$hh[1,burn], type='l', main = '', col='blue', xlab = 'Iteration') 
plot(param_phy$hh[3,burn], type='l', main = '', col='blue', xlab = 'Iteration') 


## plot(param_phy$hh[1,burn]/param_phy$hh[2,burn], type='l', main = '', col='blue', xlab = 'Iteration') 
## plot(param_phy$hh[3,burn]/param_phy$hh[4,burn], type='l', main = '', col='blue', xlab = 'Iteration') 
dev.off()


pdf(paste0(dataset, '_rho_post.pdf'))
par(mfrow=c(1,1))
hist(rowMeans(param_phy$w[,burn]), breaks=50, col='grey',main ='', ylab = 'Frequency',xlab = expression(paste('Averge posterior of ', rho)), border=TRUE)
box()
dev.off()

pdf(paste0(dataset, '_gamma_post.pdf'))
hist(rowMeans(param_phy$y[,burn]), breaks=35, col='grey',main ='', ylab = 'Frequency',xlab = expression(paste('Average posterior of ', gamma)))
box()
dev.off()

pdf(paste0(dataset, '_eta_post.pdf'))
hist(param_phy$eta[burn], breaks=30, col='grey',main ='', ylab = 'Frequency',xlab = expression(paste('Posterior of ', eta)), border=TRUE)
box()
dev.off()

pdf(paste0(dataset, '_boxplot_rho_10.pdf'))
##aux = order(rowMeans(param_phy$w[,burn]), decreasing = TRUE);aux[1:10]
aux = apply(param_phy$w[, burn],1, median )
aux = order(aux, decreasing = T);aux[1:10]

boxplot(data.frame(t(param_phy$w[aux[1:10], burn])), outline=F, col='grey', names = paste(1:10), ylab= 'Value', xlab = 'Ordered parameters' )
dev.off()


pdf(paste0(dataset, '_boxplot_gamma_10.pdf'))
##aux = order(rowMeans(param_phy$w[,burn]), decreasing = TRUE);aux[1:10]
aux = apply(param_phy$y[, burn],1, median )
aux = order(aux, decreasing = T);aux[1:10]

boxplot(data.frame(t(param_phy$y[aux[1:10], burn])), outline=F, col='grey', names = paste(1:10), ylab= 'Value', xlab = 'Ordered parameters' )
dev.off()

## Histogram
aux = getMean(param_phy)
com_pa = 1*(com>0)
P = 1-exp(-outer(aux$y, aux$w))
P = 1-exp(-outer(aux$y, aux$w)*((phy_dist^aux$eta)%*%com_pa))
pdf('affinity-single.pdf')
roc = rocCurves(Z=com_pa, Z_cross = com_pa, P=P,plot=TRUE, bins=400, all=TRUE)
dev.off()
zz = 1*(com>0)
aux= which(colSums(zz)==1)
zz[,aux]<-0
roc = rocCurves(Z=1*(com>0), Z_cross = zz, P=P,plot=TRUE, bins=400, all=FALSE)

pdf(paste0(dataset, '_hist_obs_unk.pdf'))

colass = rgb(0,0.8,0.8,0.5)
colnoass = rgb(1,0,0.4,0.4,0.5)
hist(log(P[com>0]),col=colass,main ='', ylab = 'Probability', freq=FALSE, xlab='Log of probability', ylim = c(0,0.38), xlim=c(-12,0), breaks=25)
hist(log(P[com==0]),col=colnoass,freq=FALSE, add=TRUE, breaks=25)
legend(x='topleft', legend=c('Observed associations', 'Unobserved associations'), lwd=4, col=c(colass, colnoass))

dev.off()

Z_est = 1*(roc$P>roc$threshold)
plot_Z(Z_est)
#Z_est[com_pa==1]<-1
pdf(paste0(dataset, '_degree_est_host.pdf'))
plot_degree(1*(com>0), Z_est, type='hosts')
dev.off()
pdf(paste0(dataset, '_degree_est_parasite.pdf'))
plot_degree(1*(com>0), Z_est, type='parasites')
dev.off()
pdf(paste0(dataset, '_degree_est.pdf'))
plot_degree(1*(com>0), Z_est)
dev.off()


## Paramteres and credibility interval 
aux = order(rowMeans(param_phy$w[,burn]), decreasing = TRUE);aux[1:10]
RHO = data.frame(mu = rowMeans(param_phy$w[aux[1:3], burn]), sd= apply(param_phy$w[aux[1:3], burn], 1, sd), t(apply(param_phy$w[aux[1:3], burn],1, function(r) quantile(r,probs = c(0.05, 0.95)))))
rownames(RHO)<-paste0('rho', 1:3)
aux = order(rowMeans(param_phy$y[,burn]), decreasing = TRUE);aux[1:10]
Y = data.frame(mu = rowMeans(param_phy$y[aux[1:3], burn]), sd= apply(param_phy$y[aux[1:3], burn], 1, sd), t(apply(param_phy$y[aux[1:3], burn],1, function(r) quantile(r,probs = c(0.05, 0.95)))))
rownames(Y)<-paste0('gamma', 1:3)
ETA = c(mu = mean(param_phy$eta[burn]), sd= sd(param_phy$eta[burn]))
ETA = c(ETA , as.vector(quantile(param_phy$eta[burn],probs = c(0.05, 0.95))))
TB = rbind(RHO, Y, eta =ETA)
TB

sink(paste0(dataset, '_para.txt'))
print(TB)
print(roc$auc)
TBx=data.frame(Parameter = c('$\\rho_{(1)}$','$\\rho_{(2)}$','$\\rho_{(3)}$',
                   '$\\gamma_{(1)}$','$\\gamma_{(2)}$','$\\gamma_{(3)}$', '$\\eta$'),
    Estimate = TB$mu, sd = TB$sd,
    CI = paste0('(', round(TB$X5.,2), ', ', round(TB$X95.,2) , ')'))
print(xtable(TBx), include.rownames = FALSE,,sanitize.text.function=function(x){x})
sink()

##################################################
##################################################
##################################################
## 10 fold cross validation
## ## All files
##################################################
##################################################
## ## All files
rm(list=ls())
files = grep('10foldCV-', list.dirs(), value=TRUE)
files = paste0(files,'/param.RData')
for(data in c('eid', 'gmp')){
nn.files = grep(data, files, useBytes=TRUE, value=TRUE)
nn.files = grep('Carni', nn.files, useBytes=T,invert=T, value=T)
nn.files = grep('Rod', nn.files, useBytes=T,invert=T, value=T)
nn.files = grep('Uncertain', nn.files, useBytes=T,invert=T, value=T)
nn.files = nn.files[order(nn.files)]

gres= lapply(nn.files, function(f){
    #a =regexpr(paste0('-[A-Za-z0-9]*',data,'-.*h[0-9]{2}'), f, perl=TRUE)
    a =regexpr(paste0('[A-Za-z-0-9-]*10foldCV-[A-za-z0-9-]*',data), f, perl=TRUE)
    name = substr(f, a, a + attributes(a)$match.length-1)
    name =sub('10foldCV-', '', name)
    print(name)
    load(f)
    FPR = rowMeans(sapply(res, function(r) r$FPR))
    TPR = rowMeans(sapply(res, function(r) r$TPR))
    ## if(!grepl('-NN[0-9]{2}', f)){       
    m.auc= sapply(res, function(r) r$tb$auc)
    m.thresh=sapply(res, function(r) r$tb$thres)
    m.pred=sapply(res,function(r) r$tb$pred)
    m.hold.out=sapply(res, function(r) r$tb$hold.out)
    ## }else{
    ##     m.auc= sapply(res, function(r) r$auc)
    ##     m.thresh=sapply(res, function(r) r$thres)
    ##     m.pred=sapply(res,function(r) r$pred)
    ##     m.hold.out=sapply(res, function(r) r$hold.out)
    ## }
    name = paste0(name, round(mean(m.auc),2),'-',round(mean(m.pred),2))
    TB = data.frame(m.auc, m.thresh, m.pred, m.hold.out)
    list(name=name, graph=cbind(FPR, TPR), ana=TB)
})

pdf(paste0(data, '10foldCV.pdf'))
gnames = nn.files
## gnames= c('LS-network: full model', 'Nearest-neighbour', 'LS-network: phylogeny-only','LS-network: affinity-only','LS-network: weighted-by-counts')
gcol = c('red', 'blue', 'brown', 'cyan', 'darkgreen', 'yellow')
glty = c(3,2,1,4,6,7)
#gpch = c('+', '*', 'o', 6, 7)
gpch = c(3, 8, 1, 6, 5,9)
t = 'b'
glwd=3
i= 1
plot(gres[[i]]$graph, xlab='1-specificity', ylab = 'sensitivity', type =t, col=gcol[i], main = 'ROC Curve', xlim = c(0,1), ylim = c(0,1), lty=glty[i], lwd=glwd, pch=gpch[i])
i=2
lines(gres[[i]]$graph,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
i =3
lines(gres[[i]]$graph,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
i =4
lines(gres[[i]]$graph,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
i =5
lines(gres[[i]]$graph,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
legend('bottomright', legend = gnames, col=gcol, lty=glty,lwd=2, pch = gpch)
i =6
lines(gres[[i]]$graph,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
legend('bottomright', legend = gnames, col=gcol, lty=glty,lwd=2, pch = gpch)
dev.off()

tb = t(sapply(gres, function(r) colMeans(r$ana)))
rownames(tb)<-gnames
write.csv(file=paste0(data, '-ana-10foldCV.csv'), round(tb,4))
tb
}


##################################################
##################################################
## with GMP Carnivora
rm(list=ls())
filename = 'Uncertain-10foldCV-Carnivora'
files =  paste0(grep(filename, list.dirs(), value=TRUE)[1], '/param.RData')
load(files)


## load('/home/mo/Desktop/HP-sim/Uncertain-10CV-gmp-25-02-06h58/param.RData')
G<-paramRegularG$L[1,]
## rm(paramRegular, paramRegularG)
## save.image('Uncertain-10CV-gmp-25-02-06h58/param_simple.RData')
#source('~/Dropbox/HP-Prediction/Scripts/library.R')

W = rowMeans(sapply(res, function(r) r$param$w))
Y = rowMeans(sapply(res, function(r) r$param$y))
Eta = mean(sapply(res, function(r) r$param$eta))
g  =  mean(sapply(res, function(r) r$param$L[1]))
P = 1- exp(-outer(Y,W^Eta)*((phy_dist^Eta)%*%com) )
P1 = g*P/(1-P + g*P)
P1[com>0]<-P[com>0]
P=P1

data='gmp_carnivora'

## Histograms
P1 = P
pdf(paste0(data, 'with_g_hist_obs_unk.pdf'))
colass = rgb(0,0.8,0.8,0.5)
colnoass = rgb(1,0,0.4,0.4,0.5)
hist(log(P1[com>0]),col=colass,main ='', ylab='Density', freq=FALSE, xlab='Log of probability',ylim=c(0, 0.38), xlim=c(-12,0), breaks=25)
hist(log(P1[com==0]),col=colnoass,freq=FALSE, add=TRUE, breaks=30)
legend(x='topleft', legend=c('Observed associations', 'Unobserved associations'), lwd=4, col=c(colass, colnoass))
dev.off()

P1 = PRegular
pdf(paste0(data, '_hist_obs_unk.pdf'))
colass = rgb(0,0.8,0.8,0.5)
colnoass = rgb(1,0,0.4,0.4,0.5)
hist(log(P1[com>0]),col=colass,main ='', ylab = 'Density', freq=FALSE, xlab='Log of probability', ylim = c(0,0.38), xlim=c(-12,0), breaks=25)
hist(log(P1[com==0]),col=colnoass,freq=FALSE, add=TRUE, breaks=30)
legend(x='topleft', legend=c('Observed associations', 'Unobserved associations'), lwd=4, col=c(colass, colnoass))
dev.off()

## Plotting of Z
pdf(paste0(data, 'Z-2010.pdf'))
plot_Z(com10,'parasites', 'hosts' )
dev.off()

pdf(paste0(data, '_without_gZ-2010.pdf'))
plot_Z(1*(rocRegular.all$P> rocRegular.all$threshold),'parasites', 'hosts' )
dev.off()

roctemp= rocCurves(Z=1*(com10>0), Z_cross=com, P=P, all=TRUE, plot=FALSE, bins=400)
pdf(paste0(data, '_gZ-2010.pdf'))
plot_Z(1*(roctemp$P> roctemp$threshold),'parasites', 'hosts' )
dev.off()

## Histogram of G
pdf(paste0(data, 'hist_g.pdf'),height=4)
hist(G,freq=T, xlim=c(0,0.5),col='darkblue', xlab='Posterior estimate of g for the GMP-Carnivora database', main='')
dev.off()


## Degree Distribution
#roctemp= rocCurves(Z=1*(com10>0), Z_cross=com, P=P, all=TRUE, plot=FALSE, bins=400)

Z = 1*(com10>0)
Z_est = 1*(roctemp$P>roctemp$threshold)
#Z_est[com>0]<-1
pdf(paste0(data, '_degree_est_host.pdf'))
plot_degree(Z, Z_est, type='hosts')
dev.off()
pdf(paste0(data, '_degree_est_parasite.pdf'))
plot_degree(Z, Z_est, type='parasites')
dev.off()
pdf(paste0(data, '_degree_est.pdf'))
plot_degree(Z, Z_est)
dev.off()

## Table of analysis
pdf(paste0(data, 'ROC-g.pdf'))
gnames= c('LS-network: with g', 'LS-network: without g')
gcol = c('red', 'blue')
glty = c(1,2)
gpch = c(5,2)
t = 'b'
glwd=3
i=1
plot(roctemp$roc$FPR, roctemp$roc$TPR, xlab='1-specificity', ylab = 'sensitivity', type =t, col=gcol[i], main = 'ROC Curve', xlim = c(0,1), ylim = c(0,1), lty=glty[i], lwd=glwd, pch=gpch[i])
i=2
lines(rocRegular.all$roc$FPR, rocRegular.all$roc$TPR,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
abline(a = 0, b=1,col='black',lty=2, lwd=2)
i =3
legend('bottomright', legend = gnames, col=gcol, lty=glty,lwd=2, pch=gpch)
dev.off()

## Table of analysis
## For 10fold CV with g
library(xtable)
filename = paste0(data, '_result_table.txt')
zz = 1*(roctemp$P>roctemp$threshold)
dim =dim(com)
tot.int = sum(com10>0)
zeros =sum(com10==0)
TP =sum(zz[com10>0])
FN =sum(1-zz[com10>0])
FP =sum(zz[com10==0])
TN = sum(1-zz[com10==0])
tb = matrix(c(TP, FN, FP, TN),2,2)
tb = cbind(tb, rowSums(tb))
tb = rbind(tb, colSums(tb))
colnames(tb)<-c('Positive', 'Negative', 'Total')
rownames(tb)<-c('Predicted positive', 'Predicted negative', 'Total')
print(xtable(tb,digits=0,  align='lccc', caption=paste0(data, ' with g contengency table')), file=filename)

tbg= data.frame(ana.table(com=1*(com10>0), comCross=com, roc = roctemp))


zz = 1*(rocRegular.all$P>rocRegular.all$threshold)
dim =dim(com)
tot.int = sum(com10>0)
zeros =sum(com10==0)
TP =sum(zz[com10>0])
FN =sum(1-zz[com10>0])
FP =sum(zz[com10==0])
TN = sum(1-zz[com10==0])
colnames(tb)<-c('Positive', 'Negative', 'Total')
rownames(tb)<-c('Predicted positive', 'Predicted negative', 'Total')

print(xtable(tb,digits=0,  align='lccc', caption=paste0(data, ' without g contengency table')), file=filename, append=T)

tb= data.frame(ana.table(com=1*(com10>0), comCross=com, roc = rocRegular.all))
tb = round(rbind(tbg, tb),4)
rownames(tb)=c('model with g', 'Model without g')
tb = data.frame(tb, MuG = mean(G))
tb
print(xtable(tb, digits=4), file = filename, append=T)

write.table(x=tb, file=paste0(data, '_result_table.txt'),  append = TRUE, sep='\t', eol='\n', row.names=T, col.names=T, quote=F)


##################################################
##################################################
## with GMP ALL
rm(list=ls())
filename = 'Uncertain-10foldCV-gmp'
files =  paste0(grep(filename, list.dirs(), value=TRUE)[1], '/param.RData')
load(files)

G<-paramRegularG$L[1,]
## rm(paramRegular, paramRegularG)
## save.image('Uncertain-10CV-gmp-25-02-06h58/param_simple.RData')
##source('~/Dropbox/HP-Prediction/Scripts/library.R')

W = rowMeans(sapply(res, function(r) r$param$w))
Y = rowMeans(sapply(res, function(r) r$param$y))
Eta = mean(sapply(res, function(r) r$param$eta))
g  =  mean(sapply(res, function(r) r$param$L[1]))
P = 1- exp(-outer(Y,W^Eta)*((phy_dist^Eta)%*%com) )
P1 = g*P/(1-P + g*P)
P1[com>0]<-P[com>0]
P=P1

data = 'gmp_all'
## Histograms
P1 = P
pdf(paste0(data, 'with_g_hist_obs_unk.pdf'))
colass = rgb(0,0.8,0.8,0.5)
colnoass = rgb(1,0,0.4,0.4,0.5)
hist(log(P1[com>0]),col=colass,main ='', ylab='Density', freq=FALSE, xlab='Log of probability',ylim=c(0, 0.38), xlim=c(-12,0), breaks=25)
hist(log(P1[com==0]),col=colnoass,freq=FALSE, add=TRUE, breaks=30)
legend(x='topleft', legend=c('Observed associations', 'Unobserved associations'), lwd=4, col=c(colass, colnoass))
dev.off()

P1 = PRegular
pdf(paste0(data, '_hist_obs_unk.pdf'))
colass = rgb(0,0.8,0.8,0.5)
colnoass = rgb(1,0,0.4,0.4,0.5)
hist(log(P1[com>0]),col=colass,main ='', ylab = 'Density', freq=FALSE, xlab='Log of probability', ylim = c(0,0.38), xlim=c(-12,0), breaks=25)
hist(log(P1[com==0]),col=colnoass,freq=FALSE, add=TRUE, breaks=30)
legend(x='topleft', legend=c('Observed associations', 'Unobserved associations'), lwd=4, col=c(colass, colnoass))
dev.off()

## Plotting of Z
pdf(paste0(data, 'Z-2010.pdf'))
plot_Z(com10,'parasites', 'hosts' )
dev.off()

pdf(paste0(data, '_without_gZ-2010.pdf'))
plot_Z(1*(rocRegular.all$P> rocRegular.all$threshold),'parasites', 'hosts' )
dev.off()

roctemp= rocCurves(Z=1*(com10>0), Z_cross=com, P=P, all=TRUE, plot=FALSE, bins=400)
pdf(paste0(data, '_gZ-2010.pdf'))
plot_Z(1*(P> roctemp$threshold),'parasites', 'hosts' )
dev.off()

## Histogram of G
pdf(paste0(data, 'hist_g.pdf'), height=4)
hist(G,freq=T, xlim=c(0,0.5),col='darkblue', xlab='Posterior estimate of g for the GMP database', main='')
dev.off()


## Degree Distribution
#roctemp= rocCurves(Z=1*(com10>0), Z_cross=com, P=P, all=TRUE, plot=FALSE)

Z = 1*(com10>0)
Z_est = 1*(roctemp$P>roctemp$threshold)
#Z_est[com>0]<-1
pdf(paste0(data, '_degree_est_host.pdf'))
plot_degree(Z, Z_est, type='hosts')
dev.off()
pdf(paste0(data, '_degree_est_parasite.pdf'))
plot_degree(Z, Z_est, type='parasites')
dev.off()
pdf(paste0(data, '_degree_est.pdf'))
plot_degree(Z, Z_est)
dev.off()

## Table of analysis
pdf(paste0(data, 'ROC-g.pdf'))
gnames= c('LS-network: with g', 'LS-network: without g')
gcol = c('red', 'blue')
glty = c(1,2)
gpch = c(5,2)
t = 'b'
glwd=3
i=1
plot(roctemp$roc$FPR, roctemp$roc$TPR, xlab='1-specificity', ylab = 'sensitivity', type =t, col=gcol[i], main = 'ROC Curve', xlim = c(0,1), ylim = c(0,1), lty=glty[i], lwd=glwd, pch=gpch[i])
i=2
lines(rocRegular.all$roc$FPR, rocRegular.all$roc$TPR,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
abline(a = 0, b=1,col='black',lty=2, lwd=2)
i =3
legend('bottomright', legend = gnames, col=gcol, lty=glty,lwd=2, pch=gpch)
dev.off()

## Table of analysis
## For 10fold CV with g

library(xtable)
filename = paste0(data, '_result_table.txt')
zz = 1*(roctemp$P>roctemp$threshold)
dim =dim(com)
tot.int = sum(com10>0)
zeros =sum(com10==0)
TP =sum(zz[com10>0])
FN =sum(1-zz[com10>0])
FP =sum(zz[com10==0])
TN = sum(1-zz[com10==0])
tb = matrix(c(TP, FN, FP, TN),2,2)
tb = cbind(tb, rowSums(tb))
tb = rbind(tb, colSums(tb))
colnames(tb)<-c('Positive', 'Negative', 'Total')
rownames(tb)<-c('Predicted positive', 'Predicted negative', 'Total')
print(xtable(tb,digits=0,  align='lccc', caption=paste0(data, ' with g contengency table')), file=filename)

tbg= data.frame(ana.table(com=1*(com10>0), comCross=com, roc = roctemp))


zz = 1*(rocRegular.all$P>rocRegular.all$threshold)
dim =dim(com)
tot.int = sum(com10>0)
zeros =sum(com10==0)
TP =sum(zz[com10>0])
FN =sum(1-zz[com10>0])
FP =sum(zz[com10==0])
TN = sum(1-zz[com10==0])
colnames(tb)<-c('Positive', 'Negative', 'Total')
rownames(tb)<-c('Predicted positive', 'Predicted negative', 'Total')

print(xtable(tb,digits=0,  align='lccc', caption=paste0(data, ' without g contengency table')), file=filename, append=T)

tb= data.frame(ana.table(com=1*(com10>0), comCross=com, roc = rocRegular.all))
tb = round(rbind(tbg, tb),4)
tb = data.frame(tb, MuG = mean(G))
rownames(tb)=c('model with g', 'Model without g')
print(xtable(tb, digits=4), file = filename, append=T)

write.table(x=tb, file=paste0(data, '_result_table.txt'),  append = TRUE, sep='\t', eol='\n', row.names=T, col.names=T, quote=F)


## Top 10
n = 10
Ptop=P
Ptop[com>0]<-0
library(reshape2)
colnames(Ptop)<-cnames
rownames(Ptop)<-rnames
Ptop = melt(Ptop)
Ptop = Ptop[order(Ptop[,'value'], decreasing = T),]
Ptop[1:n,]

print(xtable(Ptop[1:n,], digits=4, align ='llll'), file = filename, append=T, include.rownames=FALSE)

## Lowest 10
n = 20
Ptop=P
Ptop[com>0]<-1
library(reshape2)
colnames(Ptop)<-cnames
rownames(Ptop)<-rnames
Ptop = melt(Ptop)
Ptop = Ptop[order(Ptop[,'value'], decreasing = F),]
Ptop[1:n,][n:1,]
print('Lowest 10', file=filename, append=T)
print(xtable(Ptop[1:n,][n:1,], digits=4, align ='llll'), file = filename, append=T, include.rownames=FALSE)


##################################################
##################################################
## with EID Rodentia
rm(list=ls())
filename = 'Uncertain-10foldCV-Rodentia'
files =  paste0(grep(filename, list.dirs(), value=TRUE)[1], '/param.RData')
load(files)

G<-paramRegularG$L[1,]
## rm(paramRegular, paramRegularG)
## save.image('Uncertain-subset-10CV-eid-25-02-13h48/param_simple.RData')

W = rowMeans(sapply(res, function(r) r$param$w))
Y = rowMeans(sapply(res, function(r) r$param$y))
Eta = mean(sapply(res, function(r) r$param$eta))
g  =  mean(sapply(res, function(r) r$param$L[1]))
P = 1- exp(-outer(Y,W^Eta)*((phy_dist^Eta)%*%com) )
P1 = g*P/(1-P + g*P)
P1[com>0]<-P[com>0]
P=P1

data = 'eid_rodentia'
## Histograms
P1 = P
pdf(paste0(data, 'with_g_hist_obs_unk.pdf'))
colass = rgb(0,0.8,0.8,0.5)
colnoass = rgb(1,0,0.4,0.4,0.5)
hist(log(P1[com>0]),col=colass,main ='', ylab='Density', freq=FALSE, xlab='Log of probability',ylim=c(0, 0.38), xlim=c(-12,0), breaks=20)
hist(log(P1[com==0]),col=colnoass,freq=FALSE, add=TRUE, breaks=30)
legend(x='topleft', legend=c('Observed associations', 'Unobserved associations'), lwd=4, col=c(colass, colnoass))
dev.off()

P1 = PRegular
pdf(paste0(data, '_hist_obs_unk.pdf'))
colass = rgb(0,0.8,0.8,0.5)
colnoass = rgb(1,0,0.4,0.4,0.5)
hist(log(P1[com>0]),col=colass,main ='', ylab = 'Density', freq=FALSE, xlab='Log of probability', ylim = c(0,0.38), xlim=c(-12,0), breaks=20)
hist(log(P1[com==0]),col=colnoass,freq=FALSE, add=TRUE, breaks=30)
legend(x='topleft', legend=c('Observed associations', 'Unobserved associations'), lwd=4, col=c(colass, colnoass))
dev.off()

## Plotting of Z
pdf(paste0(data, 'Z-2010.pdf'))
plot_Z(com10,'parasites', 'hosts' )
dev.off()

pdf(paste0(data, '_without_gZ-2010.pdf'))
plot_Z(1*(rocRegular.all$P> rocRegular.all$threshold),'parasites', 'hosts' )
dev.off()

roctemp= rocCurves(Z=1*(com10>0), Z_cross=com, P=P, all=TRUE, plot=FALSE, bins=400)
pdf(paste0(data, '_gZ-2010.pdf'))
plot_Z(1*(P> roctemp$threshold),'parasites', 'hosts' )
dev.off()

## Histogram of G
pdf(paste0(data, 'hist_g.pdf'), height=4)
hist(G,freq=T, xlim=c(0,1),col='darkblue', xlab='Posterior estimate of g for the EID-Rodentia database', main='')
dev.off()

## Degree Distribution
#roctemp= rocCurves(Z=1*(com10>0), Z_cross=com, P=P, all=TRUE, plot=FALSE,bins=400)

Z = 1*(com10>0)
Z_est = 1*(roctemp$P>roctemp$threshold)
#Z_est[com>0]<-1
pdf(paste0(data, '_degree_est_host.pdf'))
plot_degree(Z, Z_est, type='hosts')
dev.off()
pdf(paste0(data, '_degree_est_parasite.pdf'))
plot_degree(Z, Z_est, type='parasites')
dev.off()
pdf(paste0(data, '_degree_est.pdf'))
plot_degree(Z, Z_est)
dev.off()

## Table of analysis
pdf(paste0(data, 'ROC-g.pdf'))
gnames= c('LS-network: with g', 'LS-network: without g')
gcol = c('red', 'blue')
glty = c(1,2)
gpch = c(5,6)
t = 'b'
glwd=3
i=1
plot(roctemp$roc$FPR, roctemp$roc$TPR, xlab='1-specificity', ylab = 'sensitivity', type =t, col=gcol[i], main = 'ROC Curve', xlim = c(0,1), ylim = c(0,1), lty=glty[i], lwd=glwd, pch=gpch[i])
i=2
lines(rocRegular.all$roc$FPR, rocRegular.all$roc$TPR,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
abline(a = 0, b=1,col='black',lty=2, lwd=2)
i =3
legend('bottomright', legend = gnames, col=gcol, lty=glty,lwd=2, pch=gpch)
dev.off()

## Table of analysis
## For 10fold CV with g

library(xtable)
filename = paste0(data, '_result_table.txt')
zz = 1*(roctemp$P>roctemp$threshold)
dim =dim(com)
tot.int = sum(com10>0)
zeros =sum(com10==0)
TP =sum(zz[com10>0])
FN =sum(1-zz[com10>0])
FP =sum(zz[com10==0])
TN = sum(1-zz[com10==0])
tb = matrix(c(TP, FN, FP, TN),2,2)
tb = cbind(tb, rowSums(tb))
tb = rbind(tb, colSums(tb))
colnames(tb)<-c('Positive', 'Negative', 'Total')
rownames(tb)<-c('Predicted positive', 'Predicted negative', 'Total')
print(xtable(tb,digits=0,  align='lccc', caption=paste0(data, ' with g contengency table')), file=filename)

tbg= data.frame(ana.table(com=1*(com10>0), comCross=com, roc = roctemp))

zz = 1*(rocRegular.all$P>rocRegular.all$threshold)
dim =dim(com)
tot.int = sum(com10>0)
zeros =sum(com10==0)
TP =sum(zz[com10>0])
FN =sum(1-zz[com10>0])
FP =sum(zz[com10==0])
TN = sum(1-zz[com10==0])
tb = matrix(c(TP, FN, FP, TN),2,2)
tb = cbind(tb, rowSums(tb))
tb = rbind(tb, colSums(tb))
colnames(tb)<-c('Positive', 'Negative', 'Total')
rownames(tb)<-c('Predicted positive', 'Predicted negative', 'Total')

print(xtable(tb,digits=0,  align='lccc', caption=paste0(data, ' without g contengency table')), file=filename, append=T)

tb= data.frame(ana.table(com=1*(com10>0), comCross=com, roc = rocRegular.all))
tb = round(rbind(tbg, tb),4)
rownames(tb)=c('model with g', 'Model without g')
tb = data.frame(tb, MuG = mean(G))
print(xtable(tb, digits=4), file = filename, append=T)

write.table(x=tb, file=paste0(data, '_result_table.txt'),  append = TRUE, sep='\t', eol='\n', row.names=T, col.names=T, quote=F)

##################################################
##################################################
## with EID ALL
rm(list=ls())
filename = 'Uncertain-10foldCV-eid'
files =  paste0(grep(filename, list.dirs(), value=TRUE)[1], '/param.RData')
load(files)

G<-paramRegularG$L[1,]
## rm(paramRegular, paramRegularG)
## save.image('Uncertain-all-10foldCV-eid-25-02-16h06/param_simple.RData')
#source('~/Dropbox/HP-Prediction/Scripts/library.R')

W = rowMeans(sapply(res, function(r) r$param$w))
Y = rowMeans(sapply(res, function(r) r$param$y))
Eta = mean(sapply(res, function(r) r$param$eta))
g  =  mean(sapply(res, function(r) r$param$L[1]))
P = 1- exp(-outer(Y,W^Eta)*((phy_dist^Eta)%*%com) )
P1 = g*P/(1-P + g*P)
P1[com>0]<-P[com>0]
P=P1

data = 'eid_all'

## Histograms
P1 = P
pdf(paste0(data, 'with_g_hist_obs_unk.pdf'))
colass = rgb(0,0.8,0.8,0.5)
colnoass = rgb(1,0,0.4,0.4,0.5)
hist(log(P1[com>0]),col=colass,main ='', ylab='Density', freq=FALSE, xlab='Log of probability',ylim=c(0, 0.38), xlim=c(-12,0), breaks=20)
hist(log(P1[com==0]),col=colnoass,freq=FALSE, add=TRUE, breaks=30)
legend(x='topleft', legend=c('Observed associations', 'Unobserved associations'), lwd=4, col=c(colass, colnoass))
dev.off()

P1 = PRegular
pdf(paste0(data, '_hist_obs_unk.pdf'))
colass = rgb(0,0.8,0.8,0.5)
colnoass = rgb(1,0,0.4,0.4,0.5)
hist(log(P1[com>0]),col=colass,main ='', ylab = 'Density', freq=FALSE, xlab='Log of probability', ylim = c(0,0.38), xlim=c(-12,0), breaks=20)
hist(log(P1[com==0]),col=colnoass,freq=FALSE, add=TRUE, breaks=30)
legend(x='topleft', legend=c('Observed associations', 'Unobserved associations'), lwd=4, col=c(colass, colnoass))
dev.off()

## Plotting of Z
pdf(paste0(data, 'Z-2010.pdf'))
plot_Z(com10,'parasites', 'hosts' )
dev.off()

pdf(paste0(data, '_without_gZ-2010.pdf'))
plot_Z(1*(rocRegular.all$P> rocRegular.all$threshold),'parasites', 'hosts' )
dev.off()

roctemp= rocCurves(Z=1*(com10>0), Z_cross=com, P=P, all=TRUE, plot=FALSE, bins=400)
pdf(paste0(data, '_gZ-2010.pdf'))
plot_Z(1*(P> roctemp$threshold),'parasites', 'hosts' )
dev.off()

## Histogram of G
pdf(file=paste0(data, 'hist_g.pdf'), height=4)
hist(G,freq=T,xlim=c(0,0.5),col='darkblue', xlab='Posterior estimate of g for the EID database', main='')
dev.off()


## Degree Distribution
#roctemp= rocCurves(Z=1*(com10>0), Z_cross=com, P=P, all=TRUE, plot=FALSE)

Z = 1*(com10>0)
Z_est = 1*(roctemp$P>roctemp$threshold)
Z_est[com>0]<-1
pdf(paste0(data, '_degree_est_host.pdf'))
plot_degree(Z, Z_est, type='hosts')
dev.off()
pdf(paste0(data, '_degree_est_parasite.pdf'))
plot_degree(Z, Z_est, type='parasites')
dev.off()
pdf(paste0(data, '_degree_est.pdf'))
plot_degree(Z, Z_est)
dev.off()

## Table of analysis
pdf(paste0(data, 'ROC-g.pdf'))
gnames= c('LS-network: with g', 'LS-network: without g')
gcol = c('red', 'blue')
glty = c(1,2)
gpch = c(5,2)
t = 'b'
glwd=3
i=1
plot(roctemp$roc$FPR, roctemp$roc$TPR, xlab='1-specificity', ylab = 'sensitivity', type =t, col=gcol[i], main = 'ROC Curve', xlim = c(0,1), ylim = c(0,1), lty=glty[i], lwd=glwd, pch=gpch[i])
i=2
lines(rocRegular.all$roc$FPR, rocRegular.all$roc$TPR,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
abline(a = 0, b=1,col='black',lty=2, lwd=2)
i =3
legend('bottomright', legend = gnames, col=gcol, lty=glty,lwd=2, pch=gpch)
dev.off()

## Table of analysis
## For 10fold CV with g

library(xtable)
filename = paste0(data, '_result_table.txt')
zz = 1*(roctemp$P>roctemp$threshold)
dim =dim(com)
tot.int = sum(com10>0)
zeros =sum(com10==0)
TP =sum(zz[com10>0])
FN =sum(1-zz[com10>0])
FP =sum(zz[com10==0])
TN = sum(1-zz[com10==0])
tb = matrix(c(TP, FN, FP, TN),2,2)
tb = cbind(tb, rowSums(tb))
tb = rbind(tb, colSums(tb))
colnames(tb)<-c('Positive', 'Negative', 'Total')
rownames(tb)<-c('Predicted positive', 'Predicted negative', 'Total')
print(xtable(tb,digits=0,  align='lccc', caption=paste0(data, ' with g contengency table')), file=filename)

tbg= data.frame(ana.table(com=1*(com10>0), comCross=com, roc = roctemp))


zz = 1*(rocRegular.all$P>rocRegular.all$threshold)
dim =dim(com)
tot.int = sum(com10>0)
zeros =sum(com10==0)
TP =sum(zz[com10>0])
FN =sum(1-zz[com10>0])
FP =sum(zz[com10==0])
TN = sum(1-zz[com10==0])
colnames(tb)<-c('Positive', 'Negative', 'Total')
rownames(tb)<-c('Predicted positive', 'Predicted negative', 'Total')

print(xtable(tb,digits=0,  align='lccc', caption=paste0(data, ' without g contengency table')), file=filename, append=T)

tb= data.frame(ana.table(com=1*(com10>0), comCross=com, roc = rocRegular.all))
tb = round(rbind(tbg, tb),4)
tb = data.frame(tb, MuG = mean(G))
rownames(tb)=c('model with g', 'Model without g')
print(xtable(tb, digits=4), file = filename, append=T)

write.table(x=tb, file=paste0(data, '_result_table.txt'),  append = TRUE, sep='\t', eol='\n', row.names=T, col.names=T, quote=F)

## Top 10
n = 10
Ptop=P
Ptop[com10>0]<-0
library(reshape2)
colnames(Ptop)<-cnames
rownames(Ptop)<-rnames
Ptop = melt(Ptop)
Ptop = Ptop[order(Ptop[,'value'], decreasing = T),]
Ptop[1:n,]

print(xtable(Ptop[1:n,], digits=4, align ='llll'), file = filename, append=T, include.rownames=FALSE)

## lowest 10
n = 20
Ptop=P
Ptop[com>0]<-1
library(reshape2)
colnames(Ptop)<-cnames
rownames(Ptop)<-rnames
Ptop = melt(Ptop)
Ptop = Ptop[order(Ptop[,'value'], decreasing = F),]

Ptop[1:n,][n:1,]

print(xtable(Ptop[1:n,][n:1,], digits=4, align ='llll'), file = filename, append=T, include.rownames=FALSE)


##################################################
##################################################
##################################################
### 100 simulation for GMP database
load('100Sim-5foldCV-gmp-28-02-15h46/param.RData')
#load('10foldCV-nodisteid-29-02-10h17/param.RData')
n= length(sim100)
dist=phy_dist
com_pa = 1*(com>0)

ana.table<-function(com, comCross, roc, plot=FALSE){
    com = 1*(com>0)
    comCross  = 1*(comCross>0)
    zz = 1*(roc$P>roc$threshold)
    if(plot) plot_Z(zz)
    tb= data.frame(auc = roc$auc, thresh= roc$threshold,
        tot.inter = sum(com), hold.out= sum(abs(comCross - com)[com==1]),
        pred = sum(zz[com==1 & comCross==0])/sum(abs(comCross - com)[com==1]),
        pred.all = sum(zz[com==1])/sum(com),
        zeros =sum(comCross==0),
        TP =sum(zz[com>0]),
        FN =sum(1-zz[com>0]),
        FP =sum(zz[com==0]),
        TN = sum(1-zz[com==0]),
        no.hosts = dim(com)[1],
        no.parasite = dim(com)[2])
    tb
}
pdf(paste0(data, '-100sim.pdf'))
plot(sim100[[1]]$graph, xlab='1-specifity', ylab = 'sensitivity', type ='l', col='grey', lwd=1, main = '',pch =2, lty=4)
for(i in 2:n){
    lines(sim100[[i]]$graph, type ='b', col='grey', lwd=1, main = '',pch =2, lty=4)
}
abline(a = 0, b=1,col='grey',lty=2, lwd=2)
dev.off()

FPR = sapply(sim100, function(r) r$graph[,1])
TPR = sapply(sim100, function(r) r$graph[,2])
FPRQ=t(apply(FPR, 1, function(r) quantile(r, c(0.05, 0.5, 0.95))))
TPRQ=t(apply(TPR, 1, function(r) quantile(r, c(0.05, 0.5, 0.95))))

lines(cbind(FPRQ[,1], TPRQ[,1]), type ='b', col='red', lwd=1,pch =4, lty=4)
lines(cbind(FPRQ[,2], TPRQ[,2]), type ='b', col='red', lwd=1,pch =4, lty=4)
lines(cbind(FPRQ[,3], TPRQ[,3]), type ='b', col='red', lwd=1,pch =4, lty=4)

pred = unlist(sapply(sim100, function(r) r$ana['m.pred'] ))
auc = unlist(sapply(sim100, function(r) r$ana['m.auc'] ))
tb = rbind(summary(pred), summary(auc))
rownames(tb)<-c('pred', 'auc')
write.csv(tb, file='100sim_report.csv')



##
load('paper-sim-results/param-gmp.RData')

par(mfrow=c(2,1))
plot(param_phy1$hh[1,], col='blue', type='l')
lines(param_phy2$hh[1,], col='orange')
lines(param_phy3$hh[1,],col ='yellow')
plot(param_phy1$hh[3,], col='blue', type='l')
lines(param_phy2$hh[3,], col='orange')
lines(param_phy3$hh[3,],col ='yellow')
