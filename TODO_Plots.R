# TODO Plots

## ##############################
## # P HeatMap plot

## # rm(list=ls())
## load('gmp-01-05-00h46/param.RData')
## # load('eid-01-05-00h46/param.RData')

## com_pa = 1*(com>0)
## burn = param_phy$burn_in - 20000:1
## burn = burn[burn>0]
## aux = getMean(param_phy)
## com_pa = 1*(com>0)
## P = 1-exp(-outer(Y^Eta, W^Eta)*((phy_dist^Eta)%*%com_pa))

## P_new <- P
## P_new[com_pa==1] <- 1
## reverse <- nrow(P_new) : 1
## P_new <- P_new[reverse,]
## image(t(P_new), col=colorRampPalette(c("white","black"))(400))

##################################################
##################################################
### gmp and eid all dataset.
### files of the format gmp- or eid-

rm(list=ls())
library(xtable)
library(coda)
files = grep('./(gmp|eid)-', list.dirs(), value=TRUE)
files = paste0(files,'/param.RData')
files = files[order(files, decreasing=TRUE)]
rho_col <- "blue"
gamma_col <- "red"
eta_col <- "darkgreen"
g_col<-'ivory4'

for (dataset in c('gmp', 'eid')){
    f = grep(dataset, files, value=T)[1] # unique files for gmp and eid
    if(length(f)==0) next
    load(f)
    print(f)
    source('library.R')

    com_pa = 1*(com>0)
    r = which.max(rowSums(com_pa));r
    c = which.max(colSums(com_pa));c
    burn = param_phy$burn_in - 15000:1
    burn = burn[burn>0]
    range(burn)

    
    pdf(paste0('Z_', dataset,'.pdf'))
    plot_Z(1*(com>0) , xlab = 'parasites', ylab = 'hosts')
    dev.off()

    pdf(paste0('degree_', dataset,'.pdf'))
    plot_degree(com>0 + 0, host.col=gamma_col, parasite.col = rho_col)
    dev.off()


    pdf(paste0('param_mcmc_',dataset,'.pdf'))
    par(mfrow=c(3,1))
    plot(param_phy$y[r,burn], type='l', main='', col=gamma_col, xlab='Iteration', ylab='Host')
    plot(param_phy$w[c,burn], type='l', main = '', col=rho_col, xlab = 'Iteration', ylab='Parasite')
    plot(param_phy$eta[burn], type='l', main = '', col=eta_col, xlab = 'Iteration', ylab='Scale')
    dev.off()

    pdf(paste0('ACF_param_mcmc_',dataset,'.pdf'))
    par(mfrow=c(3,1))
    acf(param_phy$y[r,burn],lag.max=100, main='Host - most active')
    text(80, 0.8, paste0('Effective sample size: ',round(effectiveSize(param_phy$y[r,burn]))))
    round(effectiveSize(param_phy$y[r,burn]))
    acf(param_phy$w[c,burn],lag.max=100, main='Parasite - most active')
    text(80, 0.8, paste0('Effective sample size: ',round(effectiveSize(param_phy$w[c,burn]))))
    acf(param_phy$eta[burn],lag.max=100, main='Scale parameter')
    text(80, 0.8, paste0('Effective sample size: ',round(effectiveSize(param_phy$eta[burn]))))
    dev.off()
        
    burn = param_phy$burn_in - 10000:1

    ## HIST: Average RHO
    pdf(paste0('rho_post_',dataset,'.pdf'))
    par(mfrow=c(1,1))
    hist(unlist(param_phy$w[,burn]),breaks=80, col=rho_col,main ='', ylab = 'Proportion', xlab ='', border=TRUE, freq=FALSE)
    abline(v=quantile(param_phy$w[,burn],probs = c(0.05, 0.95)), col="red", lty=2, lwd=2)
    dev.off()
    
    ## HIST: Average GAMMA
    pdf(paste0('gamma_post_', dataset,'.pdf'))
    if(dataset=='eid'){
        yy = log(param_phy$y[, burn])
        prob = c(0.05, 0.95)
    } else {
        yy = param_phy$y[, burn]
        prob = c(0.95)
    }
    hist(unlist(yy), breaks=80, col=gamma_col,main ='', ylab = 'Proportion',freq= FALSE,xlab = '', border=TRUE)
    abline(v=quantile(yy,probs = prob), col="red", lty=2, lwd=2)
    dev.off()

    ## HIST: Average ETA
    pdf(paste0('eta_post_',dataset,'.pdf'))
    hist(param_phy$eta[burn], breaks=35, col=eta_col,main ='', ylab = 'Proportion',
         freq=FALSE, xlab = '', border=TRUE)
    abline(v=quantile(param_phy$eta[burn],probs = c(0.05, 0.95)), col="red", lty=2, lwd=2)
    dev.off()

##############################################################################
    ## BOXPLOT: Invidual RHOs
    pdf(paste0('boxplot_rho_100_',dataset,'.pdf'))
    aux = apply(param_phy$w[, burn],1, median )
    no.par = 80
    aux = order(aux, decreasing = T);
    boxplot(data.frame(t(param_phy$w[aux[1:no.par], burn])), outline=F, 
            names = paste(1:no.par), ylab= 'Value', xlab = 'Ordered parameters',
            whiskcol=rho_col, staplecol=rho_col, col=rho_col, medcol="white")## BOXPLOT: Invidual RHOs
    ## Adding mean Post and 95% coverage interval
    abline(h=mean(unlist((param_phy$w[,burn]))), col='red',lty=2, lwd=1 )
    abline(h=quantile((param_phy$w[,burn]),probs = c(0.05, 0.95)), col="red", lty=2, lwd=1)
    dev.off()

    ## BOXPLOT: Invidual GAMMAs
    pdf(paste0('boxplot_gamma_100_',dataset,'.pdf'))
    aux = apply(param_phy$y[, burn],1, median )
    no.par = 80
    aux = order(aux, decreasing = T)
    boxplot(data.frame(t(param_phy$y[aux[1:no.par], burn])), outline=F, 
            names = paste(1:no.par), ylab= 'Value', xlab = 'Ordered parameters',
            whiskcol=gamma_col, staplecol=gamma_col, col=gamma_col, medcol="white")
    abline(h=mean(unlist((param_phy$y[,burn]))), col='red',lty=2, lwd=1 )
    abline(h=quantile((param_phy$y[,burn]),probs = c(0.05, 0.95)), col="red", lty=2, lwd=1)
    dev.off()

    ## AUC tables
    aux = getMean(param_phy)
    com_pa = 1*(com>0)

    if(SIMPLERHO)
        P = 1-exp(-outer(aux$y, aux$w)*((phy_dist^aux$eta)%*%com_pa)) else 
    P = 1-exp(-outer(aux$y, aux$w^aux$eta)*((phy_dist^aux$eta)%*%com_pa))

    ## LOG-PROBABILITIES
    pdf(paste0('hist_obs_unk_', dataset,'.pdf'))
    colass = rgb(0,0.8,0.8,0.5)
    colnoass = rgb(1,0,0.4,0.4,0.5)
    hist(log(P[com>0]),col=colass,main ='',
         ylab = 'Probability', freq=FALSE, xlab='Log of probability',
         ylim = c(0,0.38), xlim=c(-12,0), breaks=30)
    hist(log(P[com==0]),col=colnoass,freq=FALSE, add=TRUE, breaks=30)
    legend(x='topleft', legend=c('Observed associations',
                            'Unobserved associations'), lwd=4, col=c(colass, colnoass))
    dev.off()


    ## DEGREE DISTRIBUTION
    roc = rocCurves(Z=com_pa, Z_cross = com_pa, P=P,plot=FALSE, bins=400, all=TRUE)
    Z_est = 1*(roc$P>roc$threshold)
    #Z_est[com_pa==1]<-0
    pdf(paste0('degree_est_host_',dataset,'.pdf'))
    plot_degree(1*(com>0), Z_est, type='hosts', host.col = gamma_col, parasite.col = rho_col)
    dev.off()
    pdf(paste0('degree_est_parasite_', dataset, '.pdf'))
    plot_degree(1*(com>0), Z_est, type='parasites',host.col = gamma_col, parasite.col = rho_col)
    dev.off()
    pdf(paste0('degree_est_',dataset,'.pdf'))
    plot_degree(1*(com>0), Z_est,host.col = gamma_col, parasite.col = rho_col)
    dev.off()


    ## Paramteres and credibility interval
    no.of.param = 5
    aux = order(rowMeans(param_phy$w[,burn]), decreasing = TRUE)
    RHO = data.frame(mu = rowMeans(param_phy$w[aux[1:no.of.param], burn]), sd= apply(param_phy$w[aux[1:no.of.param], burn], 1, sd), t(apply(param_phy$w[aux[1:no.of.param], burn],1, function(r) quantile(r,probs = c(0.05, 0.95)))))
    rownames(RHO)<-paste0('rho', 1:no.of.param)
    rownames(RHO)<- paste('$\\rho_{(',1:no.of.param, ')}$', sep='')
    aux = order(rowMeans(param_phy$y[,burn]), decreasing = TRUE)
    Y = data.frame(mu = rowMeans(param_phy$y[aux[1:no.of.param], burn]), sd= apply(param_phy$y[aux[1:no.of.param], burn], 1, sd), t(apply(param_phy$y[aux[1:no.of.param], burn],1, function(r) quantile(r,probs = c(0.05, 0.95)))))
    rownames(Y)<-paste0('gamma', 1:no.of.param)
    rownames(Y)<-paste('$\\gamma_{(',1:no.of.param, ')}$', sep='')
    ETA = c(mu = mean(param_phy$eta[burn]), sd= sd(param_phy$eta[burn]))
    ETA = c(ETA , as.vector(quantile(param_phy$eta[burn],probs = c(0.05, 0.95))))
    TB = rbind(RHO, Y,'$\\eta$'=ETA)
    TB

    print(TB)
    print(roc$auc)
    
    TBx=data.frame(Parameter = rownames(TB), Estimate = TB$mu, sd = TB$sd,
        CI = paste0('(', round(TB$X5.,2), ', ', round(TB$X95.,2) , ')'))

    if(dataset=='gmp')
        filename='tbestGMP.tex'

    if(dataset=='eid')
        filename='tbestEID.tex'

    print(xtable(TBx),
          include.rownames = FALSE,
          include.colnames = FALSE,
          sanitize.text.function=function(x){x},
          only.contents=TRUE, file =filename)

}
## DONE
##################################################
##################################################


##################################################
##################################################
##################################################
## 10 fold cross validation
## ## All files
rm(list=ls())
source('../library.R')
library(xtable)
files = grep('10foldCV-', list.dirs(), value=TRUE)
files = paste0(files,'/param.RData')
all.data = list()
wilcoxon = list()

for(data in c('eid', 'gmp')){

    nn.files = grep(data, files, useBytes=TRUE, value=TRUE)
    nn.files = grep('Carni', nn.files, useBytes=T,invert=T, value=T)
    nn.files = grep('Rod', nn.files, useBytes=T,invert=T, value=T)
    nn.files = grep('Uncertain', nn.files, useBytes=T,invert=T, value=T)
    nn.files = nn.files[order(nn.files)]

    gres= lapply(nn.files, function(f){
        ## a =regexpr(paste0('-[A-Za-z0-9]*',data,'-.*h[0-9]{2}'), f, perl=TRUE)
        a =regexpr(paste0('[A-Za-z-0-9-]*10foldCV-[A-za-z0-9-]*',data), f, perl=TRUE)
        name = substr(f, a, a + attributes(a)$match.length-1)
        #name =sub('10foldCV-', '', name)
        print(name)
        load(f)
        FPR = rowMeans(sapply(res, function(r) r$FPR))
        TPR = rowMeans(sapply(res, function(r) r$TPR))

        m.auc= sapply(res, function(r) r$tb$auc)
        m.thresh=sapply(res, function(r) r$tb$thres)
        m.pred=sapply(res,function(r) r$tb$pred)
        m.hold.out=sapply(res, function(r) r$tb$hold.out)

        P = NULL
        legend.name =NULL
        if(grepl('DistOnly', name)){
            Eta = mean(sapply(res, function(r) r$param))
            P = 1-exp(-(phy_dist^Eta) %*%(1*(com>0)))
            legend.name ='LS-network: phylogeny-only'
             
        }

        if(grepl('Weighted', name)){
            W = rowMeans(sapply(res, function(r) r$param$w))
            Y = rowMeans(sapply(res, function(r) r$param$y))
            Eta = mean(sapply(res, function(r) r$param$eta))
            
            ## if(data =='gmp'){
            ##     com[com>2]<-log(1+com)[com>2]
            ## } 
            if(data =='eid'){
                com=log(com+1)/2
            }
            U = (phy_dist^Eta) %*% com
            if(SIMPLERHO){
                P = 1-  exp(-outer(Y, W)*U)
            }else{
                P = 1-  exp(-outer(Y, W^Eta)*U)
            }
            legend.name = 'LS-network: weighted-by-counts'

        }

        if(grepl('nodist', name)){
            W = rowMeans(sapply(res, function(r) r$param$w))
            Y = rowMeans(sapply(res, function(r) r$param$y))
            P = 1-  exp(-outer(Y, W))
            legend.name = 'LS-network: affinity-only'

        }

        if(grepl('NN', name)){
            P<-NULL
            legend.name  ='Nearest-neighbour' 
        }


        if(name=='10foldCV-gmp' | name =='10foldCV-eid'){
            W = rowMeans(sapply(res, function(r) r$param$w))
            Y = rowMeans(sapply(res, function(r) r$param$y))
            Eta = mean(sapply(res, function(r) r$param$eta))

            U = (phy_dist^Eta)%*%(1*(com>0))
            
            if(SIMPLERHO){
                P = 1-  exp(-outer(Y, W)*U)
            }else{
                P = 1-  exp(-outer(Y, W^Eta)*U)
            }
            legend.name = 'LS-network: full model'
        }
        
        TB = data.frame(m.auc, m.thresh, m.pred, m.hold.out)
        list(name=name, legend.name = legend.name,graph=cbind(FPR, TPR), ana=TB, P=P)
    })



    ord.model = c('LS-network: full model','LS-network: affinity-only','LS-network: phylogeny-only','LS-network: weighted-by-counts','Nearest-neighbour')

    pdf(paste0('10foldCV-',data,'.pdf'))
    gnames = sapply(gres, function(r) r[['legend.name']])
    ord.plot = sapply(ord.model, function(r) which(r==gnames))
    #gnames = nn.files
    ## gnames= c('LS-network: full model', 'Nearest-neighbour', 'LS-network: phylogeny-only','LS-network: affinity-only','LS-network: weighted-by-counts')
    gcol = c('black', 'lemonchiffon4', 'brown', 'cyan', 'darkgreen', 'yellow')
    glty = c(3,2,1,4,6,7)
    ## gpch = c('+', '*', 'o', 6, 7)
    gpch = c(3, 8, 1, 6, 5,9)
    t = 'b'
    glwd=3
    i= 1
    plot(gres[[i]]$graph, xlab='1-specificity', ylab = 'sensitivity', type =t, col=gcol[i], main = 'ROC Curve', xlim = c(0,1), ylim = c(0,1), lty=glty[i], lwd=glwd, pch=gpch[i])
    for (i in 2:length(gres)){
        lines(gres[[i]]$graph,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
    }
    legend('bottomright', legend = gnames[ord.plot], col=gcol[ord.plot], lty=glty[ord.plot],lwd=2, pch = gpch[ord.plot])
    dev.off()

    ### Plotting the interaction posterior matrix 
    for( i in 1:length(gres)){
        if(!is.null(gres[[i]]$P)){
            pdf(paste0('Z_', gres[[i]]$name, '.pdf'))
            plot_Z(1*(gres[[i]]$P > mean(gres[[i]]$ana$m.thresh) ),xlab = 'parasites', ylab = 'hosts')
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
}

xx = cbind(all.data[['gmp']][,c(1,2,3)], all.data[['eid']][,c(2,3)])
rownames(xx)<-NULL
xx = xx[sapply(ord.model, function(r) which(r==xx[,'names'])),]

print(xtable(xx),
      include.rownames = FALSE,
      include.colnames = FALSE,
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

rm(list=ls())
filenames = paste0('Uncertain-10foldCV-',c('gmp','Carnivora','eid','Rodentia'))
all.data = list()
rho_col <- "blue"
gamma_col <- "red"
eta_col <- "darkgreen"
g_col<-'ivory4'

wil=list()
wilcoxon.local.test<-function(x,y){
    d = x-y
    r = rank(abs(d))
    T = min(sum(r[d>0]) + 0.5*sum(r[d==0]), sum(r[d<0])+0.5*sum(r[d==0]))
    N = length(r)
    z = (T - 0.25*N*(N+1))/(sqrt((1/24)*N*(N+1)*(2*N+1)))
    2*pnorm(-abs(z))
}


for(data in filenames){
    print(data)
    aux = grep(data, list.dirs(), value=TRUE)
    print(aux)
    if(length(aux)==0) next
    load(paste0(aux[1], '/param.RData'))
    source('../library.R')
    com=1*(com>0)
    ## Creating probability matrix

    aux = sapply(res, function(r) r[['withG']]['param'])
    G = mean(sapply(aux, function(r) r[['g']]))
    W = rowMeans(sapply(aux, function(r) r[['w']]))
    Y = rowMeans(sapply(aux, function(r) r[['y']]))
    Eta = mean(sapply(aux, function(r) r[['eta']]))

    if(SIMPLERHO){
        P = 1- exp(-outer(Y,W)*((phy_dist^Eta)%*%com) )
    }else{
        P = 1- exp(-outer(Y,W^Eta)*((phy_dist^Eta)%*%com) )
    }
       
    Pg = G*P/(1-P + G*P)
    Pg[com>0]<-P[com>0]
    
    aux = sapply(res, function(r) r[['withOutG']]['param'])
    W = rowMeans(sapply(aux, function(r) r[['w']]))
    Y = rowMeans(sapply(aux, function(r) r[['y']]))
    Eta = mean(sapply(aux, function(r) r[['eta']]))

    if(SIMPLERHO){
        Pnog = 1- exp(-outer(Y,W)*((phy_dist^Eta)%*%com) )
    }else{
        Pnog = 1- exp(-outer(Y,W^Eta)*((phy_dist^Eta)%*%com) )
    }
    ## AUC for Wilcoxon
    wil[[data]]<-cbind(all = wilcoxon.local.test(unlist(sapply(res, function(r)
        r[['withG']][['tb.all']]['auc'])),
                           unlist(sapply(res, function(r) r[['withOutG']][['tb.all']]['auc']))),
                       subset = wilcoxon.local.test(unlist(sapply(res, function(r)
                           r[['withG']][['tb']]['auc'])),
                           unlist(sapply(res, function(r) r[['withOutG']][['tb']]['auc']))))

    ## Histograms
    P1 = Pg
    pdf(paste0('with_g_hist_obs_unk_', data,'.pdf'))
    colass = rgb(0,0.8,0.8,0.5)
    colnoass = rgb(1,0,0.4,0.4,0.5)
    hist(log(P1[com10>0]),col=colass,main ='', ylab='Density', freq=FALSE, xlab='Log of probability',ylim=c(0, 0.38), xlim=c(-15,0), breaks=25)
    hist(log(P1[com10==0]),col=colnoass,freq=FALSE, add=TRUE, breaks=30)
    legend(x='topleft', legend=c('Observed associations', 'Unobserved associations'), lwd=4, col=c(colass, colnoass))
    dev.off()
    
    P1 = Pnog
    pdf(paste0('hist_obs_unk_', data,'.pdf'))
    colass = rgb(0,0.8,0.8,0.5)
    colnoass = rgb(1,0,0.4,0.4,0.5)
    hist(log(P1[com10>0]),col=colass,main ='', ylab = 'Density', freq=FALSE, xlab='Log of probability', ylim = c(0,0.38), xlim=c(-15,0), breaks=25)
    hist(log(P1[com10==0]),col=colnoass,freq=FALSE, add=TRUE, breaks=30)
    legend(x='topleft', legend=c('Observed associations', 'Unobserved associations'), lwd=4, col=c(colass, colnoass))
    dev.off()
    
    ## Plotting of Z
    pdf(paste0('Z-2010_', data,'.pdf'))
    plot_Z(com10,'parasites', 'hosts' )
    dev.off()

    rocNoG = rocCurves(Z=1*(com10>0), Z_cross=com, P=Pnog, all=TRUE, plot=FALSE, bins=400)
    pdf(paste0('without_gZ-2010_',data,'.pdf'))
    plot_Z(1*(rocNoG$P> rocNoG$threshold),'parasites', 'hosts' )
    dev.off()
    
    rocG= rocCurves(Z=1*(com10>0), Z_cross=com, P=Pg, all=TRUE, plot=FALSE, bins=400)
    pdf(paste0('gZ-2010_',data,'.pdf'))
    plot_Z(1*(rocG$P> rocG$threshold),'parasites', 'hosts' )
    dev.off()
    
    ## Histogram of G
    pdf(paste0('hist_g_',data,'.pdf'),height=4)
    hist(paramRegularG$g,freq=T,col=g_col, xlab='Posterior estimate of g for the GMP-Carnivora database', main='', breaks=40, xlim=c(0,0.6))
    abline(v=quantile(paramRegularG$g,probs = c(0.05, 0.95)), col="red", lty=2, lwd=2)
    dev.off()
    
    
    ## Degree Distribution
    Z = 1*(com10>0)
    Z_est = 1*(rocG$P>rocG$threshold)
    pdf(paste0('degree_est_host_', data,'.pdf'))
    plot_degree(Z, Z_est, type='hosts',host.col = gamma_col, parasite.col = rho_col)
    dev.off()
    pdf(paste0('degree_est_parasite_', data,'.pdf'))
    plot_degree(Z, Z_est, type='parasites',host.col = gamma_col, parasite.col = rho_col)
    dev.off()
    pdf(paste0('degree_est_',data,'.pdf'))
    plot_degree(Z, Z_est,host.col = gamma_col, parasite.col = rho_col)
    dev.off()
    
    ## Table of analysis
    pdf(paste0('ROC-g_', data,'.pdf'))
    gnames= c('LS-network: with g', 'LS-network: without g')
    gcol = c('black', 'ivory4')
    glty = c(1,2)
    gpch = c(5,2)
    t = 'b'
    glwd=3
    i=1
    plot(rocG$roc$FPR, rocG$roc$TPR, xlab='1-specificity', ylab = 'sensitivity', type =t, col=gcol[i], main = 'ROC Curve', xlim = c(0,1), ylim = c(0,1), lty=glty[i], lwd=glwd, pch=gpch[i])
    i=2
    lines(rocNoG$roc$FPR, rocNoG$roc$TPR,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
    abline(a = 0, b=1,col='black',lty=2, lwd=2)
    legend('bottomright', legend = gnames, col=gcol, lty=glty,lwd=2, pch=gpch)
    dev.off()
    
    ## Table of analysis
    ## For 10fold CV with g
    library(xtable)
    filename = paste0('result_table_',data,'.txt')
    
    tbg= data.frame(ana.table(com=1*(com10>0), comCross=com, roc = rocG))
    tb= data.frame(ana.table(com=1*(com10>0), comCross=com, roc = rocNoG))
    tb = round(rbind(tbg, tb),3)
    rownames(tb)=c('Model with $g$', 'Model without $g$')
    tb = data.frame(tb, MuG = mean(G))
    tb
  
    TBx = tb[,c('auc','pred','pred.all', 'MuG')]
    TBx = data.frame(TBx)
    
    all.data[[sub('Uncertain-10foldCV-', '', data)]]<-TBx
   
}

## AUC table: 
AUCtb = sapply(all.data, function(r) round(r[['auc']]/100,3))
AUCtb = cbind(Model=c('with $g$', 'without $g$'), AUCtb)

print(xtable(AUCtb, digits=3),
      include.colnames = TRUE,
      include.rownames = FALSE,
      sanitize.text.function=function(x){x},
      only.contents=TRUE,file='tb-GMP-EID-AUC.tex')
    
## Prediction from Subsets and all sets
## Prediction table with subet
PRED =  sapply(all.data, function(r) paste0('(',round(r[['pred']],3),') ', round(r[['pred.all']],3)) )
colnames(PRED)<-names(all.data)
PRED = cbind(Model=c('with $g$', 'without $g$'), PRED)
rownames(AUCtb)=c('with $g$', 'without $g$')
print(xtable(PRED),
      include.colnames = TRUE,
      include.rownames = FALSE,
      sanitize.text.function=function(x){x},
      only.contents=TRUE, file='tb-GMP-EID-Pred.tex')



MuG = sapply(all.data, function(r) round(r[['MuG']],4))[1,]
MuG = data.frame(t(MuG))
print(xtable(MuG, digits = 3),
      include.colnames = TRUE,
      include.rownames = TRUE,
      sanitize.text.function=function(x){x},
      only.contents=TRUE,file='tb-GMP-EID-MuG.tex')

WIL = data.frame(do.call('rbind', wil))
rownames(WIL)<-sub('Uncertain-10foldCV-', '', names(wil))

print(xtable(WIL, digits = 3),
      include.colnames = TRUE,
      include.rownames = TRUE,
      sanitize.text.function=function(x){x},
      only.contents=TRUE,file='tbWilcoxon-Uncer.tex')

### DONE
##################################################
