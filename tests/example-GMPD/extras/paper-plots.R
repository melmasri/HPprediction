s##################################################
##################################################
### gmp and eid all dataset.
### files of the format gmp- or eid-

rm(list=ls())
library(xtable)
library(coda)
files = grep('(ALL)-', list.dirs(recurse  =FALSE), value=TRUE)
files = paste0(files,'/param.RData')
files = files[order(files, decreasing=TRUE)]
rho_col <- "blue"
gamma_col <- "red"
eta_col <- "darkgreen"
g_col<-'ivory4'

dataset='gmp'

f = grep(dataset, files, value=T, ignore.case= TRUE)[1] # unique files for gmp and eid
if(length(f)==0) next
load(f)
print(f)
source('../../library.R')

com_pa = 1*(com>0)
r = which.max(rowSums(com_pa));r
c = which.max(colSums(com_pa));c
burn = param$burn_in - 15000:1
burn = burn[burn>0]
range(burn)
##com = com[order(host.order),]
##param$y = param$y[order(host.order),]

pdf(paste0('Z_', dataset,'.pdf'))
plot_Z(1*(com[order(host.order),]>0) , xlab = 'parasites', ylab = 'hosts')
dev.off()

pdf(paste0('degree_', dataset,'.pdf'))
plot_degree(com>0 + 0, host.col=gamma_col, parasite.col = rho_col)
dev.off()


pdf(paste0('param_mcmc_',dataset,'.pdf'))

par(mfrow=c(3,1))
plot(param$y[r,burn], type='l', main='', col=gamma_col, xlab='Iteration', ylab='Host',
     cex.lab=1.5)
plot(param$w[c,burn], type='l', main = '', col=rho_col, xlab = 'Iteration', ylab='Parasite', cex.lab=1.5)
plot(param$eta[burn], type='l', main = '', col=eta_col, xlab = 'Iteration', ylab='Tree scaling parameter', cex.lab=1.5)

dev.off()

pdf(paste0('ACF_param_mcmc_',dataset,'.pdf'))

par(mfrow=c(3,1))
acf(param$y[r,burn],lag.max=100, main='Host - most active')
text(80, 0.8, paste0('Effective sample size: ',round(effectiveSize(param$y[r,burn]))))
round(effectiveSize(param$y[r,burn]))
acf(param$w[c,burn],lag.max=100, main='Parasite - most active')
text(80, 0.8, paste0('Effective sample size: ',round(effectiveSize(param$w[c,burn]))))
round(effectiveSize(param$w[c,burn]))
acf(param$eta[burn],lag.max=100, main='Scale parameter')
text(80, 0.8, paste0('Effective sample size: ',round(effectiveSize(param$eta[burn]))))
round(effectiveSize(param$eta[burn]))

dev.off()
        

## HIST: Average RHO
pdf(paste0('rho_post_',dataset,'.pdf'))

par(mfrow=c(1,1))
hist(log(unlist(param$w[,burn])),breaks=80,
     col=rho_col,main ='', ylab = 'Proportion', xlab ='', border=TRUE, freq=FALSE)
abline(v=quantile(log(param$w[,burn]),probs = c(0.05, 0.95)), col="red", lty=2, lwd=2)

dev.off()
    
## HIST: Average GAMMA
pdf(paste0('gamma_post_', dataset,'.pdf'))

hist(log(param$y[, burn]), breaks=80, col=gamma_col,main ='', ylab = 'Proportion',freq= FALSE,xlab = '', border=TRUE)
abline(v=quantile(log(param$y[, burn]),probs = c(0.05,0.95)), col="red", lty=2, lwd=2)

dev.off()

    ## HIST: Average ETA
pdf(paste0('eta_post_',dataset,'.pdf'))
hist(param$eta[burn], breaks=35, col=eta_col,main ='', ylab = 'Proportion',
     freq=FALSE, xlab = '', border=TRUE)
abline(v=quantile(param$eta[burn],probs = c(0.05, 0.95)), col="red", lty=2, lwd=2)
dev.off()

##############################################################################
## BOXPLOT: Invidual RHOs
pdf(paste0('boxplot_rho_100_',dataset,'.pdf'))

aux = apply(param$w[, burn],1, median )
no.par = 80
aux = order(aux, decreasing = T);
boxplot(data.frame(t(param$w[aux[1:no.par], burn])), outline=F, 
        names = paste(1:no.par), ylab= 'Value', xlab = 'Ordered parameters',
        whiskcol=rho_col, staplecol=rho_col, col=rho_col, medcol="white")## BOXPLOT: Invidual RHOs
## Adding mean Post and 95% coverage interval
abline(h=mean(unlist((param$w[,burn]))), col='red',lty=2, lwd=1 )
abline(h=quantile((param$w[,burn]),probs = c(0.05, 0.95)), col="red", lty=2, lwd=1)

dev.off()

## BOXPLOT: Invidual GAMMAs
pdf(paste0('boxplot_gamma_100_',dataset,'.pdf'))

aux = apply(param$y[, burn],1, median )
no.par = 80
aux = order(aux, decreasing = T)
boxplot(data.frame(t(param$y[aux[1:no.par], burn])), outline=F, 
        names = paste(1:no.par), ylab= 'Value', xlab = 'Ordered parameters',
        whiskcol=gamma_col, staplecol=gamma_col, col=gamma_col, medcol="white")
abline(h=mean(unlist((param$y[,burn]))), col='red',lty=2, lwd=1 )
abline(h=quantile((param$y[,burn]),probs = c(0.05, 0.95)), col="red", lty=2, lwd=1)

dev.off()

    ## AUC tables
aux = getMean(param)
com_pa = 1*(com>0)
library(ape)
library(geiger)
dd = 1/cophenetic(rescale(tree, 'EB', aux$eta))
diag(dd)<-0
Z= unname(com[host.order,])
Z = 1*(Z>0)
pdist = dd %*% Z
P = 1-exp(-outer(aux$y, aux$w)*pdist)

## LOG-PROBABILITIES
pdf(paste0('hist_obs_unk_', dataset,'.pdf'))
colass = rgb(0,0.8,0.8,0.5)
colnoass = rgb(1,0,0.4,0.4,0.5)
hist(log(P[Z==0]),col=colnoass,main ='',
     ylab = 'Probability', freq=FALSE, xlab='Log of probability', ylim=c(0, 0.3))
hist(log(P[Z>0]),col=colass,freq=FALSE, add=TRUE)
legend(x='topright', legend=c('Observed associations',
                         'Unobserved associations'), lwd=4, col=c(colass, colnoass))
dev.off()
    
## Paramteres and credibility interval
no.of.param = 5
aux = order(rowMeans(param$w[,burn]), decreasing = TRUE)
RHO = data.frame(mu = rowMeans(param$w[aux[1:no.of.param], burn]), sd= apply(param$w[aux[1:no.of.param], burn], 1, sd), t(apply(param$w[aux[1:no.of.param], burn],1, function(r) quantile(r,probs = c(0.05, 0.95)))))
rownames(RHO)<-paste0('rho', 1:no.of.param)
rownames(RHO)<- paste('$\\rho_{(',1:no.of.param, ')}$', sep='')
aux = order(rowMeans(param$y[,burn]), decreasing = TRUE)
Y = data.frame(mu = rowMeans(param$y[aux[1:no.of.param], burn]), sd= apply(param$y[aux[1:no.of.param], burn], 1, sd), t(apply(param$y[aux[1:no.of.param], burn],1, function(r) quantile(r,probs = c(0.05, 0.95)))))
rownames(Y)<-paste0('gamma', 1:no.of.param)
rownames(Y)<-paste('$\\gamma_{(',1:no.of.param, ')}$', sep='')
ETA = c(mu = mean(param$eta[burn]), sd= sd(param$eta[burn]))
ETA = c(ETA , as.vector(quantile(param$eta[burn],probs = c(0.05, 0.95))))
TB = rbind(RHO, Y,'$\\eta$'=ETA)
TB

print(TB)
##print(roc$auc)

TBx=data.frame(Parameter = rownames(TB), Estimate = TB$mu, sd = TB$sd,
    CI = paste0('(', round(TB$X5.,2), ', ', round(TB$X95.,2) , ')'))
colnames(TBx)<-c('Parameter', 'Estimate', 'standard dev', '95\\% credible interval')
if(dataset=='gmp')
    filename='tbestGMP.tex'

if(dataset=='eid')
    filename='tbestEID.tex'

print(xtable(TBx),
      include.rownames = FALSE,
      include.colnames = TRUE,
      sanitize.text.function=function(x){x},
      only.contents=TRUE, file =filename)

roc= rocCurves(Z , Z, TRUE, P, 400,TRUE)
## DONE
##################################################
##################################################


##################################################
##################################################
##################################################
## 10 fold cross validation
## ## All files
rm(list=ls())
source('~/Github/HP-ICM/extras/library.R')
library(xtable)
library(ape)
library(geiger)

files = grep('GMPD', list.dirs(), value=TRUE)
files = paste0(files,'/param.RData')
all.data = list()
wilcoxon = list()

data = 'gmp'

nn.files = grep(data, files, useBytes=TRUE, value=TRUE, ignore.case= TRUE)
nn.files = grep('Carni', nn.files, useBytes=T,invert=T, value=T)
nn.files = grep('Rod', nn.files, useBytes=T,invert=T, value=T)
nn.files = grep('Uncertain', nn.files, useBytes=T,invert=T, value=T)
nn.files = nn.files[order(nn.files)]

gres= lapply(nn.files, function(f){
    ## a =regexpr(paste0('-[A-Za-z0-9]*',data,'-.*h[0-9]{2}'), f, perl=TRUE)
    a =regexpr('[A-Za-z-0-9-]*com[A-z-.]*', f, perl=TRUE)
    name = substr(f, a, a + attributes(a)$match.length-1)
    ##name =sub('10foldCV-', '', name)
    print(name)
    load(f)
    FPR = rowMeans(sapply(res, function(r) r$FPR))
    TPR = rowMeans(sapply(res, function(r) r$TPR))
    indices= lof(com, TRUE)
    m.auc = sapply(res, function(r) r$tb$auc)/100
    m.thresh = sapply(res, function(r) r$tb$thres)
    m.pred = sapply(res,function(r) r$tb$pred)*100
    m.hold.out = sapply(res, function(r) r$tb$hold.out)
    
    P = NULL
    legend.name =NULL
    if(grepl('DistOnly', name, ignore.case = TRUE)){
        Eta = mean(sapply(res, function(r) r$param$eta))
        pdist = 1/cophenetic(rescale(tree, 'EB', Eta))
        diag(pdist)<-0
        pdist = pdist %*% com_pa
        pdist[pdist==0]<-10e10
        P = 1- exp(-pdist)
        P = P[,indices]
        legend.name ='LS-net: phylogeny-only'
    }
    
    if(grepl('Weighted', name, ignore.case = TRUE)){
        W = rowMeans(sapply(res, function(r) r$param$w))
        Y = rowMeans(sapply(res, function(r) r$param$y))
        Eta = mean(sapply(res, function(r) r$param$eta))
        pdist = 1/cophenetic(rescale(tree, 'EB', Eta))
        diag(pdist)<-0
        pdist = pdist %*% (log(com+1)/2)
        P = 1-  exp(-outer(Y, W)*pdist)
        P = P[,indices]
        legend.name = 'LS-net: weighted-by-counts'
    }

    if(grepl('NN', name)){
        legend.name = 'Nearest-neighbour'
        P = NULL
    }
    
    if(grepl('affinity', name, ignore.case=TRUE)){
        W = rowMeans(sapply(res, function(r) r$param$w))
        Y = rowMeans(sapply(res, function(r) r$param$y))
        P = 1-  exp(-outer(Y, W))
        P = P[, indices]
        legend.name = 'LS-net: affinity-only'
        }
    
    if(grepl('fold', name, ignore.case=TRUE)){
        W = rowMeans(sapply(res, function(r) r$param$w))
        Y = rowMeans(sapply(res, function(r) r$param$y))
        Eta = mean(sapply(res, function(r) r$param$eta))

        pdist = 1/cophenetic(rescale(tree, 'EB', Eta))
        diag(pdist)<-0
        pdist = pdist %*% com_pa
        pdist[pdist==0]<-1
        P = 1-  exp(-outer(Y, W)*pdist)
        P = P[,indices]
        legend.name = 'LS-net: full model'
    }
    TB = data.frame(m.auc, m.thresh, m.pred, m.hold.out)
    list(name=name, legend.name = legend.name,graph=cbind(FPR, TPR), ana=TB, P=P)
})

ord.model = c('LS-net: full model','LS-net: affinity-only','LS-net: phylogeny-only','LS-net: weighted-by-counts','Nearest-neighbour')

pdf(paste0('10foldCV-',data,'.pdf'))
gnames = sapply(gres, function(r) r[['legend.name']])
ord.plot = sapply(ord.model, function(r) which(r==gnames))
ord.plot = unlist(ord.plot)
gcol = c('black', 'lemonchiffon4', 'brown', 'cyan', 'darkgreen', 'yellow')
glty = c(3,2,1,4,6,7)
## gpch = c('+', '*', 'o', 6, 7)
gpch = c(3, 8, 1, 6, 5,9)
t = 'b'
glwd=3
i= 1
plot(gres[[i]]$graph, xlab='1-specificity', ylab = 'sensitivity', type =t, col=gcol[i], main = 'ROC Curve', xlim = c(0,1), ylim = c(0,1),
     lty=glty[i], lwd=glwd,  pch=gpch[i],
     cex.lab=1.5)
if(length(gres)>1)
    for (i in 2:length(gres)){
        lines(gres[[i]]$graph,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
    }
legend('bottomright', legend = gnames[ord.plot], col=gcol[ord.plot], lty=glty[ord.plot],lwd=2, pch = gpch[ord.plot], pt.cex=1, cex=1.5)
dev.off()

    ### Plotting the interaction posterior matrix 
for( i in 1:length(gres)){
    if(!is.null(gres[[i]]$P)){
        pdf(paste0('Z_', gres[[i]]$name, '.pdf'))
        plot_Z(1*(gres[[i]]$P > mean(gres[[i]]$ana$m.thresh) ),xlab = '', ylab = '')
        dev.off()
    }
}
    
tb = t(sapply(gres, function(r) colMeans(r$ana)))
rownames(tb)<-gnames
write.csv(file=paste0('Ana-10foldCV-',data,'.csv'), round(tb,4))
    
## LATEX TABLE

names = sub('/param.RData', '', gnames)
names = sub('./', '', names)
names = sub('10foldCV-', '', names)
    
TBx = data.frame(names, tb[,'m.auc'], tb[,'m.pred'])
all.data[[data]]<-TBx

tb = t(sapply(gres, function(r) r$ana[,'m.auc']))
rownames(tb)<-gnames
wilcoxon[[data]]<-tb


xx = cbind(all.data[['gmp']][,c(1,2,3)])
rownames(xx)<-NULL
xx = xx[unlist(sapply(ord.model, function(r) which(r==xx[,'names']))),]
colnames(xx)<-c('Model', 'AUC', 'Prediction')

print(xtable(xx, digits=3),
      include.rownames = FALSE,
      include.colnames = TRUE,
      sanitize.text.function=function(x){x},
      only.contents=TRUE, file='tbAUC.tex')



wilcoxon.local.test<-function(x,y){
    d = x-y
    r = rank(abs(d))
    T = min(sum(r[d>0]) + 0.5*sum(r[d==0]), sum(r[d<0])+0.5*sum(r[d==0]))
    N = length(r)
    z = (T - 0.25*N*(N+1))/(sqrt((1/24)*N*(N+1)*(2*N+1)))
    2*pnorm(-abs(z))
}

ll=list()

for(n in names(wilcoxon)){
    dd = wilcoxon[[n]]
    grid = expand.grid(1:nrow(dd),1:nrow(dd))
    p = cbind(grid, p = apply(grid,1, function(r)
        wilcoxon.local.test(dd[r[1],],dd[r[2],])))

    wil = matrix(0, nrow(dd), nrow(dd))
    wil[cbind(p[,1], p[,2])]<-round(p[,3],5)
    wil[upper.tri(wil)]<-0
    rownames(wil)<-sub('LS-network: ', '', rownames(dd))
    ll[[n]]<-wil
}

print(xtable(cbind(ll[['gmp']], ll[['eid']]),digits=3),
      include.rownames = TRUE,
      include.colnames = FALSE,
      sanitize.text.function=function(x){x},
      only.contents=TRUE,
      file='tbWilcoxon-10fold-AUC.tex' )

## DONE
##################################################
##################################################


##################################################
##################################################
## Uncertainty
##################################################
setwd('../')

rm(list=ls())
files = grep('uncertain', list.dirs(recurse=FALSE), value=TRUE, ignore.case = TRUE)
all.data = list()
rho_col <- "blue"
gamma_col <- "red"
eta_col <- "darkgreen"
g_col<-'ivory4'
all.top.m = list()
wil=list()

wilcoxon.local.test<-function(x,y){
    d = x-y
    r = rank(abs(d))
    T = min(sum(r[d>0]) + 0.5*sum(r[d==0]), sum(r[d<0])+0.5*sum(r[d==0]))
    N = length(r)
    z = (T - 0.25*N*(N+1))/(sqrt((1/24)*N*(N+1)*(2*N+1)))
    2*pnorm(-abs(z))
}
library(ape)
library(geiger)

for(f  in files){

    print(f)
    a =regexpr('./UNCERTAIN-[A-za-z0-9-]*-year', f, perl=TRUE)
    name = substr(f, a, a + attributes(a)$match.length-1)
    name = gsub('./UNCERTAIN-', '', name)
    name = gsub('-year', '', name)
        
    if(length(f)==0) next

    setwd(f)
    load('param.RData')
    ##load(paste0(f[1], '/param.RData'))
    com=1*(com>0)
    ## Creating probability matrix
    
    aux = sapply(res, function(r) r[['withG']]['param'])
    G = mean(sapply(aux, function(r) r[['g']]))
    W = rowMeans(sapply(aux, function(r) r[['w']]))
    Y = rowMeans(sapply(aux, function(r) r[['y']]))
    Eta = mean(sapply(aux, function(r) r[['eta']]))

    pdist = 1/cophenetic(rescale(tree, 'EB', Eta))
    diag(pdist)<-0
    pdist = pdist %*% com
    pdist[pdist==0] <- 1
    P = 1-  exp(-outer(Y, W)*pdist)
    
    Pg = G*P/(1-P + G*P)
    Pg[com>0]<-P[com>0]
    
    aux = sapply(res, function(r) r[['withOutG']]['param'])
    W = rowMeans(sapply(aux, function(r) r[['w']]))
    Y = rowMeans(sapply(aux, function(r) r[['y']]))
    Eta = mean(sapply(aux, function(r) r[['eta']]))

    pdist = 1/cophenetic(rescale(tree, 'EB', Eta))
    diag(pdist)<-0
    pdist = pdist %*% com
    pdist[pdist==0] <-1
    Pnog = 1-  exp(-outer(Y, W)*pdist)


    
    ## AUC for Wilcoxon
    wil[[name]]<-cbind(all = wilcoxon.local.test(unlist(sapply(res, function(r)
        r[['withG']][['tb.all']]['auc'])),
                           unlist(sapply(res, function(r) r[['withOutG']][['tb.all']]['auc']))),
                       subset = wilcoxon.local.test(unlist(sapply(res, function(r)
                           r[['withG']][['tb']]['auc'])),
                           unlist(sapply(res, function(r) r[['withOutG']][['tb']]['auc']))))

    ## Histograms
    P1 = Pg
    pdf(paste0('with_g_hist_obs_unk_', name,'.pdf'))
    colass = rgb(0,0.8,0.8,0.5)
    colnoass = rgb(1,0,0.4,0.4,0.5)
    hist(log(P1[com10>0]),col=colass,main ='', ylab='Density', freq=FALSE, xlab='Log of probability',ylim=c(0, 0.38), xlim=c(-11,0), breaks=30)
    hist(log(P1[com10==0]),col=colnoass,freq=FALSE, add=TRUE)
    legend(x='topleft', legend=c('Observed associations', 'Unobserved associations'), lwd=4, col=c(colass, colnoass),pt.cex=1, cex=1.5) 
    dev.off()
    
    P1 = Pnog
    pdf(paste0('hist_obs_unk_', name,'.pdf'))
    colass = rgb(0,0.8,0.8,0.5)
    colnoass = rgb(1,0,0.4,0.4,0.5)
    hist(log(P1[com10>0]),col=colass,main ='', ylab = 'Density', freq=FALSE, xlab='Log of probability', ylim = c(0,0.38), xlim=c(-10,0), breaks=30)
    hist(log(P1[com10==0]),col=colnoass,freq=FALSE, add=TRUE)
    legend(x='topleft', legend=c('Observed associations', 'Unobserved associations'), lwd=4, col=c(colass, colnoass),pt.cex=1, cex=1.5)
    dev.off()
    
    ind = lof(com10, TRUE)
    ## Plotting of Z
    pdf(paste0('Z-2010_', name,'.pdf'))
    plot_Z(com10[,ind],'parasites', 'hosts' )
    dev.off()

    rocNoG = rocCurves(Z=1*(com10>0), Z_cross=com, P=Pnog, all=TRUE, plot=FALSE, bins=400)
    pdf(paste0('without_gZ-2010_',name,'.pdf'))
    plot_Z((1*(rocNoG$P> rocNoG$threshold))[,ind],'parasites', 'hosts' )
    dev.off()
    
    rocG= rocCurves(Z=1*(com10>0), Z_cross=com, P=Pg, all=TRUE, plot=FALSE, bins=400)
    pdf(paste0('gZ-2010_',name,'.pdf'))
    plot_Z((1*(rocG$P> rocG$threshold))[,ind],'', '' )
    dev.off()
    
    ## Histogram of G
    pdf(paste0('hist_g_',name,'.pdf'),height=4)
    hist(paramG$g,freq=F,col=g_col, xlab='Posterior estimate of g', main='', breaks=100, xlim=c(0,0.4), cex.lab=1.5)
    abline(v=quantile(paramG$g,probs = c(0.05, 0.95)), col="red", lty=2, lwd=2)
    dev.off()

    #acf(c(paramG$g),lag.max=100, main='Scale parameter')
    #text(80, 0.8, paste0('Effective sample size: ',round(effectiveSize(paramG$g))))
    
    cbind(Model=c('with G', 'without G'), rbind( round(ana.table(com10, com, rocG), 4), round(ana.table(com10, com, rocNoG), 4)))

    ## Degree Distribution G
    Z = 1*(com10>0)
    Z_est = 1*(rocG$P>rocG$threshold)
    pdf(paste0('degree_est_host_G_', name,'.pdf'))
    plot_degree(Z, Z_est, type='hosts',host.col = gamma_col, parasite.col = rho_col)
    dev.off()
    pdf(paste0('degree_est_parasite_G_', name,'.pdf'))
    plot_degree(Z, Z_est, type='parasites',host.col = gamma_col, parasite.col = rho_col)
    dev.off()
    pdf(paste0('degree_est_G_',name,'.pdf'))
    plot_degree(Z, Z_est,host.col = gamma_col, parasite.col = rho_col)
    dev.off()
    
    
    ## Degree Distribution No G
    Z = 1*(com10>0)
    Z_est = 1*(rocNoG$P>rocNoG$threshold)
    pdf(paste0('degree_est_host_NoG_', name,'.pdf'))
    plot_degree(Z, Z_est, type='hosts',host.col = gamma_col, parasite.col = rho_col)
    dev.off()
    pdf(paste0('degree_est_parasite_NoG_', name,'.pdf'))
    plot_degree(Z, Z_est, type='parasites',host.col = gamma_col, parasite.col = rho_col)
    dev.off()
    pdf(paste0('degree_est_NoG_',name,'.pdf'))
    plot_degree(Z, Z_est,host.col = gamma_col, parasite.col = rho_col)
    dev.off()
    
    ## Table of analysis
    pdf(paste0('ROC-g_', name,'.pdf'))
    gnames= c('LS-net full model: with g', 'LS-net full model: without g')
    gcol = c('black', 'ivory4')
    glty = c(1,2)
    gpch = c(5,2)
    t = 'b'
    glwd=3
    i=1
    plot(rocG$roc$FPR, rocG$roc$TPR, xlab='1-specificity', ylab = 'sensitivity', type =t, col=gcol[i], main = 'ROC Curve', xlim = c(0,1), ylim = c(0,1), lty=glty[i], lwd=glwd, pch=gpch[i], cex.lab=1.5)
    i=2
    lines(rocNoG$roc$FPR, rocNoG$roc$TPR,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
    abline(a = 0, b=1,col='black',lty=2, lwd=2)
    legend('bottomright', legend = gnames, col=gcol, lty=glty,lwd=2, pch=gpch,
           pt.cex=1, cex=1.5)
    dev.off()
    
    ## Table of analysis
    ## For 10fold CV with g
    library(xtable)
    filename = paste0('result_table_',name,'.txt')
    
    tbg= data.frame(ana.table(com=1*(com10>0), comCross=com, roc = rocG))
    tb= data.frame(ana.table(com=1*(com10>0), comCross=com, roc = rocNoG))
    tb = rbind(tbg, tb)
    rownames(tb)=c('Model with $g$', 'Model without $g$')
    tb = data.frame(tb, MuG = mean(G))
    tb
  
    TBx = tb[,c('auc','pred','pred.all', 'MuG')]
    TBx = data.frame(TBx)
    
    all.data[[sub('Uncertain-10foldCV-', '', name)]]<-TBx

    ## top m
    m = 4000
    Z = 1*(com10>0)
    ## Full with G
    ord.pg = order(rocG$P, decreasing=TRUE)
    Zpg = 1*(rocG$P>rocG$threshold)
    ## Full without G
    ord.p = order(rocNoG$P, decreasing=TRUE)
    Zp = 1*(rocNoG$P>rocNoG$threshold)

    topm = t(sapply(1:m, function(r) c(actual = sum(Z[ord.pg[1:r]]),
        withG= sum(Z[ord.pg[1:r]]*Zpg[ord.pg[1:r]]),
        withOutG= sum(Z[ord.p[1:r]]*Zp[ord.p[1:r]])))) 


    topm = data.frame(topm)
    all.top.m[[sub('Uncertain-10foldCV-', '', name)]]<-rbind(topm[100,],topm[1000,])
    
    pdf(paste0('TopM_', name,'.pdf'))

    plot(x=1:m,y = topm[,'withG'], xlab='Number of validated pairwise interactions', ylab = 'Number of recovered pairwise interactions',  col='red',lty=1, type='l', lwd=2,
         cex.lab=1.5)
    lines(x=1:m,y = topm[,'withOutG'], lty=3, type='l', lwd=2, col='red')
    gnames = c('Full model with uncertainty', 'Full model')
    gcol = c('red','red')
    glwd = c(2, 2)
    glty = c(1, 3)

    
    lines(x=1:m, y=1:m, lty=5, lwd=1, col='black')
    gnames = c(gnames, 'x=y')
    gcol = c(gcol, 'black')
    glwd = c(glwd,1)
    glty = c(glty,5)
    legend('bottomright', legend = gnames, col=gcol, lty=glty,lwd=glwd,
           pt.cex=1,cex=1.5)
    dev.off()


    


clean.names<-function(aux){
    aux= sub('Carnivora-gmp','gmp-Carnivora', aux, ignore.case=TRUE)
    aux= sub('Rodentia-eid','eid-Rodentia', aux, ignore.case=TRUE)
    aux= sub('gmp','GMPD', aux, ignore.case=TRUE)
    aux= sub('eid','EID2', aux, ignore.case=TRUE)
    aux = sub('Uncertain-10foldCV-', '',aux)
    aux = sub('-affin', ' Affinity', aux, ignore.case=TRUE)
    aux = sub('-dist', ' Phylogeny', aux, ignore.case=TRUE)
    aux
}

## AUC table: 
AUCtb = sapply(all.data, function(r) round(r[['auc']]/100,3))
AUCtb = cbind(Model=c('with $g$', 'without $g$'), AUCtb)
aux = colnames(AUCtb)
colnames(AUCtb)<-clean.names(colnames(AUCtb))

gmp.only = grepl('(model|gmp)', colnames(AUCtb), ignore.case=TRUE)
print(xtable(AUCtb[,gmp.only], digits=3),
      include.colnames = TRUE,
      include.rownames = FALSE,
      sanitize.text.function=function(x){x},
      only.contents=TRUE,file='tb-GMP-AUC.tex')

major.only = !grepl('(Aff|Phy)', colnames(AUCtb), ignore.case=TRUE)
print(xtable(AUCtb[,major.only], digits=3),
      include.colnames = TRUE,
      include.rownames = FALSE,
      sanitize.text.function=function(x){x},
      only.contents=TRUE,file='tb-GMP-EID-AUC.tex')
    
## Prediction from Subsets and all sets
## Prediction table with subet
PRED =  sapply(all.data, function(r) paste0('(',round(100*r[['pred']],3),') ', round(100*r[['pred.all']],3)) )
colnames(PRED)<-clean.names(colnames(PRED))
PRED = cbind(Model=c('with $g$', 'without $g$'), PRED)

gmp.only = grepl('(model|gmp)', colnames(PRED), ignore.case=TRUE)
major.only = !grepl('(Aff|Phy)', colnames(PRED), ignore.case=TRUE)

print(xtable(PRED[,major.only]),
      include.colnames = TRUE,
      include.rownames = FALSE,
      sanitize.text.function=function(x){x},
      only.contents=TRUE, file='tb-GMP-EID-Pred.tex')

print(xtable(PRED[,gmp.only]),
      include.colnames = TRUE,
      include.rownames = FALSE,
      sanitize.text.function=function(x){x},
      only.contents=TRUE, file='tb-GMP-Pred.tex')


MuG = sapply(all.data, function(r) round(r[['MuG']],4))[1,]
names(MuG)<-clean.names(names(MuG))
MuG = data.frame(t(MuG))
print(xtable(MuG, digits = 3),
      include.colnames = TRUE,
      include.rownames = TRUE,
      sanitize.text.function=function(x){x},
      only.contents=TRUE,file='tb-GMP-EID-MuG.tex')

WIL = data.frame(do.call('rbind', wil))
rownames(WIL)<-clean.names(names(wil))

print(xtable(WIL, digits = 3),
      include.colnames = TRUE,
      include.rownames = TRUE,
      sanitize.text.function=function(x){x},
      only.contents=TRUE,file='tbWilcoxon-Uncer.tex')

    setwd('../')
}
## TOPM

### DONE
##################################################


## Plot Tree
pdf('tree-GMPD.pdf')
plot(tree, font=1, cex=0.5,label.offset=1)
dev.off()



