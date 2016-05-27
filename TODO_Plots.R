# TODO Plots

rm(list=ls())
load('param.RData')

rho_col <- "blue"
gamma_col <- "red"
eta_col <- "darkgreen"


# Original Data
com_pa = 1*(com>0)

burn = param_phy$burn_in - 20000:1
burn = burn[burn>0]

# pdf(paste0(dataset, '_Z.pdf'))
plot_Z(1*(com>0) , xlab = 'parasites', ylab = 'hosts')
# dev.off()

# HIST: Average RHO
# pdf(paste0(dataset, '_rho_post.pdf'))
hist(rowMeans(param_phy$w[,burn]), breaks=50, col=rho_col,main ='', ylab = 'Frequency',
	xlab = expression(paste('Averge posterior of ', rho)), border=TRUE, xlim=c(23,32))
abline(v=quantile(rowMeans(param_phy$w[,burn]),probs = c(0.05, 0.95)), col="red", lty=2, lwd=2)
# box()
# dev.off()

# HIST: Average GAMMA
# pdf(paste0(dataset, '_gamma_post.pdf'))
hist(rowMeans(param_phy$y[,burn]), breaks=35, col=gamma_col,main ='', ylab = 'Frequency',
	xlab = expression(paste('Averge posterior of ', gamma)), border=TRUE)
abline(v=quantile(rowMeans(param_phy$y[,burn]),probs = 0.95), col="red", lty=2, lwd=2)
# box()
# dev.off()

# HIST: Average ETA
# pdf(paste0(dataset, '_eta_post.pdf'))
hist(param_phy$eta[burn], breaks=30, col=eta_col,main ='', ylab = 'Frequency',
	xlab = expression(paste('Posterior of ', eta)), border=TRUE, xlim=c(1.32, 1.48))
abline(v=quantile(param_phy$eta[burn],probs = c(0.05, 0.95)), col="red", lty=2, lwd=2)
# box()
# dev.off()


# TODO
# HIST: g for both datasets


#######################################
# BOXPLOT: Invidual RHOs
# pdf(paste0(dataset, '_boxplot_rho_100.pdf'))
aux = apply(param_phy$w[, burn],1, median )
aux = order(aux, decreasing = T);aux[1:100]
boxplot(data.frame(t(param_phy$w[aux[1:100], burn])), outline=F, 
	names = paste(1:100), ylab= 'Value', xlab = 'Ordered parameters',
	whiskcol=rho_col, staplecol=rho_col, col=rho_col, medcol="white")# BOXPLOT: Invidual RHOs
# dev.off()


# BOXPLOT: Invidual GAMMAs
# pdf(paste0(dataset, '_boxplot_gamma_100.pdf'))
aux = apply(param_phy$y[, burn],1, median )
aux = order(aux, decreasing = T);aux[1:100]
boxplot(data.frame(t(param_phy$y[aux[1:100], burn])), outline=F, 
	names = paste(1:100), ylab= 'Value', xlab = 'Ordered parameters',
	whiskcol=gamma_col, staplecol=gamma_col, col=gamma_col, medcol="white")
# dev.off()



##############################
# P HeatMap plot

# rm(list=ls())
load('gmp-01-05-00h46/param.RData')
# load('eid-01-05-00h46/param.RData')

com_pa = 1*(com>0)
burn = param_phy$burn_in - 20000:1
burn = burn[burn>0]
aux = getMean(param_phy)
com_pa = 1*(com>0)
P = 1-exp(-outer(Y^Eta, W^Eta)*((phy_dist^Eta)%*%com_pa))

P_new <- P
P_new[com_pa==1] <- 1
reverse <- nrow(P_new) : 1
P_new <- P_new[reverse,]
image(t(P_new), col=colorRampPalette(c("white","black"))(400))
image(t(P_new), col=colorRampPalette(c("white","blue","brown"))(400))




######################################
# Tables
require(xtable)
rm(list=ls())
load("/home/max/GitHub/HP-prediction/April28/Uncertain-10foldCV-gmp-27-04-20h54/param.RData")
ls()

filename <- "test.tex"
roctemp= rocCurves(Z=1*(com10>0), Z_cross=com, P=P, all=TRUE, plot=FALSE, bins=400)
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
print(xtable(tb,digits=0,  align='lccc'), file=filename)



##################################################
##################################################
### gmp and eid all dataset.
### files of the format gmp- or eid-

rm(list=ls())
files = grep('./(gmp|eid)-', list.dirs(), value=TRUE)
files = paste0(files,'/param.RData')
files = files[order(files, decreasing=TRUE)]



rho_col <- "blue"
gamma_col <- "red"
eta_col <- "darkgreen"


for (dataset in c('gmp', 'eid')){
    f = grep(dataset, files, value=T)[1] # unique files for gmp and eid
    if(length(f)==0) next

    load(f)

    com_pa = 1*(com>0)
    r = which.max(rowSums(com_pa));r
    c = which.max(colSums(com_pa));c
    burn = param_phy$burn_in - 15000:1
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
    plot(param_phy$y[r,burn], type='l', main='', col=rho_col, xlab='Iteration', ylab='Host parameter')
    plot(param_phy$w[c,burn], type='l', main = '', col=gamma_col, xlab = 'Iteration', ylab='Parasite parameter')
    plot(param_phy$eta[burn], type='l', main = '', col=eta_col, xlab = 'Iteration', ylab='Scaling parameter')
    dev.off()

    burn = param_phy$burn_in - 10000:1


    ## HIST: Average RHO
    pdf(paste0(dataset, '_rho_post.pdf'))
    par(mfrow=c(1,1))
    hist(rowMeans(param_phy$w[,burn]),breaks=80, col=rho_col,main ='', ylab = 'Frequency', xlab = expression(paste('Averge posterior of ', rho)), border=TRUE)
    abline(v=quantile(rowMeans(param_phy$w[,burn]),probs = c(0.05, 0.95)), col="red", lty=2, lwd=2)
    dev.off()
    
    ## HIST: Average GAMMA
    pdf(paste0(dataset, '_gamma_post.pdf'))
    hist(rowMeans(param_phy$y[,burn]), breaks=30, col=gamma_col,main ='', ylab = 'Frequency',
         xlab = expression(paste('Averge posterior of ', gamma)), border=TRUE)
    abline(v=quantile(rowMeans(param_phy$y[,burn]),probs = 0.95), col="red", lty=2, lwd=2)
    dev.off()

    ## HIST: Average ETA
    pdf(paste0(dataset, '_eta_post.pdf'))
    hist(param_phy$eta[burn], breaks=30, col=eta_col,main ='', ylab = 'Frequency',
         xlab = expression(paste('Posterior of ', eta)), border=TRUE, xlim=c(1.32, 1.48))
    abline(v=quantile(param_phy$eta[burn],probs = c(0.05, 0.95)), col="red", lty=2, lwd=2)
    dev.off()

##############################################################################
    ## BOXPLOT: Invidual RHOs
    pdf(paste0(dataset, '_boxplot_rho_100.pdf'))
    aux = apply(param_phy$w[, burn],1, median )
    aux = order(aux, decreasing = T);aux[1:100]
    boxplot(data.frame(t(param_phy$w[aux[1:100], burn])), outline=F, 
            names = paste(1:100), ylab= 'Value', xlab = 'Ordered parameters',
            whiskcol=rho_col, staplecol=rho_col, col=rho_col, medcol="white")## BOXPLOT: Invidual RHOs
    dev.off()

    ## BOXPLOT: Invidual GAMMAs
    pdf(paste0(dataset, '_boxplot_gamma_100.pdf'))
    aux = apply(param_phy$y[, burn],1, median )
    aux = order(aux, decreasing = T);aux[1:100]
    boxplot(data.frame(t(param_phy$y[aux[1:100], burn])), outline=F, 
            names = paste(1:100), ylab= 'Value', xlab = 'Ordered parameters',
            whiskcol=gamma_col, staplecol=gamma_col, col=gamma_col, medcol="white")
    dev.off()

    ## AUC tables
    aux = getMean(param_phy)
    com_pa = 1*(com>0)

    if(SIMPLERHO)
        P = 1-exp(-outer(Y, W)*((phy_dist^Eta)%*%com_pa)) else 
    P = 1-exp(-outer(Y, W^Eta)*((phy_dist^Eta)%*%com_pa))

    ## LOG-PROBABILITIES
    pdf(paste0(dataset, '_hist_obs_unk.pdf'))
    colass = rgb(0,0.8,0.8,0.5)
    colnoass = rgb(1,0,0.4,0.4,0.5)
    hist(log(P[com>0]),col=colass,main ='',
         ylab = 'Probability', freq=FALSE, xlab='Log of probability',
         ylim = c(0,0.38), xlim=c(-12,0), breaks=25)
    hist(log(P[com==0]),col=colnoass,freq=FALSE, add=TRUE, breaks=25)
    legend(x='topleft', legend=c('Observed associations',
                            'Unobserved associations'), lwd=4, col=c(colass, colnoass))
    dev.off()


    ## DEGREE DISTRIBUTION
    roc = rocCurves(Z=com_pa, Z_cross = com_pa, P=P,plot=TRUE, bins=400, all=TRUE)
    Z_est = 1*(roc$P>roc$threshold)
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
source('library.R')
files = grep('10foldCV-', list.dirs(), value=TRUE)
files = paste0(files,'/param.RData')

all.data = list()

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
        name =sub('10foldCV-', '', name)
        print(name)
        load(f)
        FPR = rowMeans(sapply(res, function(r) r$FPR))
        TPR = rowMeans(sapply(res, function(r) r$TPR))

        m.auc= sapply(res, function(r) r$tb$auc)
        m.thresh=sapply(res, function(r) r$tb$thres)
        m.pred=sapply(res,function(r) r$tb$pred)
        m.hold.out=sapply(res, function(r) r$tb$hold.out)

        P = NULL
        if(!is.null(res[[1]]$param)){
            W = rowMeans(sapply(res, function(r) r$param$w))
            Y = rowMeans(sapply(res, function(r) r$param$y))
            if(!is.null(res[[1]]$param$eta)){
                Eta = mean(sapply(res, function(r) r$param$eta))
            } else Eta=1
            
            SIMPLERHO = TRUE
            dist = phy_dist
            
            if(grepl('nodist', name)){
                U = 1
            }else  if(grepl('weighted', name)){
                if(data =='gmp'){
                    com[com>2]<-log(1+com)[com>2]
                } 
                if(data =='eid'){
                    com=log(com+1)/2
                }
                U = (dist^Eta) %*% com
            }else U = (dist^Eta)%*%(1*(com>0))

            if(!grepl('NN', name)){
                if(SIMPLERHO){
                    P = 1-  exp(-outer(Y, W)*U)
                }else{
                    P = 1-  exp(-outer(Y, W))
                }
            }
        }
        
        TB = data.frame(m.auc, m.thresh, m.pred, m.hold.out)
        list(name=name, graph=cbind(FPR, TPR), ana=TB, P=P)
    })


    pdf(paste0(data, '10foldCV.pdf'))
    gnames = nn.files
    ## gnames= c('LS-network: full model', 'Nearest-neighbour', 'LS-network: phylogeny-only','LS-network: affinity-only','LS-network: weighted-by-counts')
    gcol = c('red', 'blue', 'brown', 'cyan', 'darkgreen', 'yellow')
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
    legend('bottomright', legend = gnames, col=gcol, lty=glty,lwd=2, pch = gpch)
    dev.off()

    ### Plotting the interaction posterior matrix 
    for( i in 1:length(gres)){
        if(!is.null(gres[[i]]$P)){
            pdf(paste0('Z_', gres[[i]]$name, '.pdf'))
            plot_Z(1*(gres[[i]]$P > mean(gres[[i]]$ana$m.thresh) ))
            dev.off()
        }
    }
    
    tb = t(sapply(gres, function(r) colMeans(r$ana)))
    rownames(tb)<-gnames
    write.csv(file=paste0(data, '-ana-10foldCV.csv'), round(tb,4))
    
    ## LATEX TABLE

    names = sub('/param.RData', '', gnames)
    names = sub('./', '', names)
    names = sub('10foldCV-', '', names)
    
    TBx = data.frame(names, tb[,'m.auc'], tb[,'m.pred'])
    print(xtable(TBx),
          include.rownames = FALSE,
          include.colnames = FALSE,
          sanitize.text.function=function(x){x},
          only.contents=TRUE, file=paste0('tbAUC', data, '.tex'))
}

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
for(data in filenames){
    aux = grep(data, list.dirs(), value=TRUE)
    if(length(aux)==0) next
    load(paste0(aux[1], '/param.RData'))

    com=1*(com>0)
    ## Creating probability matrix
    G<-paramRegularG$L[1,]
    W = rowMeans(sapply(res, function(r) r$param$w))
    Y = rowMeans(sapply(res, function(r) r$param$y))
    Eta = mean(sapply(res, function(r) r$param$eta))
    g  =  mean(sapply(res, function(r) r$param$L[1]))

    if(SIMPLERHO){
        P = 1- exp(-outer(Y,W)*((phy_dist^Eta)%*%com) )
    }else{
        P = 1- exp(-outer(Y,W^Eta)*((phy_dist^Eta)%*%com) )
    }

    P1 = g*P/(1-P + g*P)
    P1[com>0]<-P[com>0]
    P=P1
    
    
    

    ## Histograms
    P1 = P
    pdf(paste0(data, 'with_g_hist_obs_unk.pdf'))
    colass = rgb(0,0.8,0.8,0.5)
    colnoass = rgb(1,0,0.4,0.4,0.5)
    hist(log(P1[com>0]),col=colass,main ='', ylab='Density', freq=FALSE, xlab='Log of probability',ylim=c(0, 0.38), xlim=c(-15,0), breaks=25)
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
    Z = 1*(com10>0)
    Z_est = 1*(roctemp$P>roctemp$threshold)
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

    legend('bottomright', legend = gnames, col=gcol, lty=glty,lwd=2, pch=gpch)
    dev.off()
    
    ## Table of analysis
    ## For 10fold CV with g
    library(xtable)
    filename = paste0(data, '_result_table.txt')
    
    tb= data.frame(ana.table(com=1*(com10>0), comCross=com, roc = rocRegular.all))
    tb = round(rbind(tbg, tb),4)
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

### DONE
##################################################
