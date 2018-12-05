### ==================================================
### Analysis functions
### ==================================================


### ==================================================
### Main functions
### ==================================================
plot_Z<-function(Z, xlab, ylab, ...){
    ## ploting interaction matrix as a binary image
    if(missing(ylab)) ylab = 'hosts'
	if(missing(xlab)) xlab = 'parasites'
    par(mar = c(5,5,1,1)+0.1)
	image(1:ncol(Z), 1:nrow(Z), t(Z[nrow(Z):1,]),
		col = c('white', 'black'), ylab=ylab, xlab=xlab,
          useRaster=TRUE,srt=45, axes=FALSE,cex.lab=3)
           

    a = 100*round(ncol(Z)*0.25 /100,0)
	axis(1, at = a*0:(ceiling(ncol(Z)/a)), cex.axis= 2)
    b = 100*round(nrow(Z)*0.25 /100,0)
	axis(2, at = b*0:(ceiling(nrow(Z)/b)),cex.axis = 2)
}

lof<-function(Z, indices = FALSE){
    ## Given a binary matrix Z. Where the rows is fixed
    ## A function that left orders the matrix sequentially from row 1 to n
    ## based on first appearance of columns.
    if(min(range(Z))<0) stop('Range is less that 0.')

    active_col <- apply(Z,1,function(r) which(as.vector(r)>0))
	bank = active_col[[1]]
	for(i in 1:nrow(Z)){
			a = setdiff(active_col[[i]], bank)
			bank=c(bank,a)
		}
    if(indices) bank else  Z[,bank]
}

rocCurves<-function(Z.test,Z.train,P, plot=TRUE,bins=400, all=FALSE){
    ## Computes the ROC curves as x = FPR and y = TPR
    ##  FPR = FP/(FP +TN)
    ##  TPR = TP/(TP+FN)
    ## Arguments:
    ## plot = TRUE to actually plot the ROC curve
    ## bin  = the number of intervals (0,1) is divided into
    ## all = TRUE, when TPR and FPR are calculated on the whole dataset, FALSE on the held out portion only
    ## Returns:
    ## auc = the maximum AUC value
    ## threshold = the threshold where (P> threshold) has maximum AUC
    ## roc = a matrix of (threshold, FPR, TPR)
    Z.test= 1*(Z.test>0)
    Z.train = 1*(Z.train>0)
    
    u = seq(0, 1, length.out = bins)    # binning
    Z1 = (Z.test - Z.train)==1
    if(all) Z1 = Z.test==1
    m = sum(1*(Z1))
    Z2 = Z.test==0
    n = sum(1*Z2)
    aux = sapply(u, function(r){
                     aux = 1*(P>=r)
                     TP = sum(aux[Z1])            # True positive
                     FN = m - TP                  # False negative
                     FP = sum(aux[Z2])           # False positive
                     TN = n - FP                # True negative
                     c(TP=TP, FN=FN, FP=FP, TN=TN)
                 })
    FPR = aux['FP',]/(aux['FP',] + aux['TN',])
    TPR = aux['TP',]/(aux['TP',] + aux['FN',])
    roc = data.frame(u=u, FPR = FPR, TPR = TPR)
    max.point =   which.min(abs(roc$TPR-1) + roc$FPR)
    threshold = roc$u[max.point]
    n = length(u)
    auc =  0.5*t(abs(FPR[2:n] - FPR[2:n -1]))%*% (TPR[2:n] + TPR[2:n -1])
    auc = round(100*auc,2)
    if(plot){
        plot(FPR, TPR, xlab='FPR (1-specifity)', ylab = 'TPR (sensitivity)',
             type ='l', col='red', lwd=2, main = paste('AUC ', auc),
             xlim = c(0,1), ylim = c(0,1), pch =6, lty=4)
        abline(a = 0, b=1,col='black',lty=2, lwd=2)
    }
    list(auc = auc, threshold = threshold,roc=roc)
}

topPairs<-function(P,Z,topX=20){
    ## Returning pairs with highest posterior probability
    require(reshape2)
    P[Z>0]<--1
    aux =   melt(P)
    aux = aux[order(aux$value,decreasing=TRUE),]
    colnames(aux)<-c('Hosts', 'Parasites', 'p')
    aux[1:topX,]
}


cross.validate.fold<-function(Z, n= 10, min.per.col = 1, missing.pattern=c('random','prop.to.col.sums')){
    ## n-fold cross validation
    ## Returns a matrix of 3 columns, the first two are the (row,col) index of the pair,
    ## the third is the group
    missing.pattern = tolower(missing.pattern[1])
    if(max(range(Z))>1) Z[Z>0]<-1
    pairs = which(Z==1, arr.ind=T)
    colnames(pairs)<-c('row', 'col')
    
    if(length(which(colSums(Z)<min.per.col))>0){
        aux = which(pairs[,'col'] %in% which(colSums(Z)<min.per.col))
        if(length(aux))
            pairs = pairs[-aux,]
    }
    
    colm = pmax(colSums(Z) -min.per.col , 0)
    size = floor(sum(colm)/n)
    gr = rep(size, n)
    if(sum(colm) %% size!=0)
        gr[n] =  gr[n] + sum(colm) %% size

    pair.list = list()
    for(i in 1:n){
        bank=c()
        for(k in 1:gr[i]){
            a = which(colm>0)
            if(missing.pattern=='random')
                b = a[sample(length(a),1)] else
            if (missing.pattern=='prop.to.col.sums')
                b = a[sample(length(a),1, prob=colm[a]/sum(colm[a]))] else
            stop('missing pattern has to be specified from selection!')
            bank = c(bank, b)
            colm[b] = colm[b]-1
        }
        pair.list[[i]]<-bank
    }

    
    gr.list= list()
    bank= c()
    for(i in 1:n){
        a= table(pair.list[[i]])
        gr.rows = unlist(sapply(1:length(a), function(r){
            b = which(pairs[,'col']== as.numeric(names(a[r])))
            b =setdiff(b, bank)
            b[sample.int(length(b), a[r])]
        }))
        bank = c(bank, gr.rows)
        gr.list[[i]]<-cbind(gr.rows, i)
    }

    aux = do.call('rbind', gr.list)
    pairs = cbind(pairs[aux[,1], ],gr= aux[,2])
    
    print(sprintf("Actual cross-validation rate is %0.3f" , table(pairs[,'gr'])/sum(1*(Z>0))))
    pairs[order(pairs[,'gr']),]
    
}

cross.validate.set<-function(Z, rate= 0.2){
    ##  Retunrs a Z with a percentage of one's turned 0
    ## only cells with rows of more than one interaction and columns with more than 2
    if(max(range(Z))>1) Z[Z>0]<-1
    Zo = Z
    pairs = which(Z==1, arr.ind=T)
    n = floor(rate*nrow(pairs))
    sampled = rep(0, nrow(pairs))
    repeat{
        r = sample(which(sampled==0),1)
        if(sum(Z[pairs[r,1],])>1 & sum(Z[,pairs[r,2]])>2){
            Z[pairs[r,1], pairs[r,2]]<-0
            sampled[r]<-1
            n = n -1
        }
        if(n==0) break
    }
    print(sprintf("Actual cross-validation rate is %0.3f" ,sum(Zo-Z)/sum(Zo)))
    Z
}

### ==================================================
### Secondary functions
### ==================================================

ana.table<-function(Z, ZCross, P, roc, plot=FALSE){
    Z = 1*(Z>0)
    ZCross  = 1*(ZCross>0)
    Zpost = 1*(P>roc$threshold)
    data.frame(auc = roc$auc/100, thresh= roc$threshold,
               tot.ones = sum(Z), held.out.ones= sum(abs(ZCross - Z)[Z==1]),
               pred.held.out.ones = 100*sum(Zpost[Z==1 & ZCross==0])/sum(abs(ZCross - Z)[Z==1]),
               pred.tot.ones = sum(Zpost[Z==1])/sum(Z)*100)
    
}

plot_degree <- function(Z, Z_est, type='both', host.col='blue', parasite.col='red'){
    ## Plots the degree distribution per marginal on a bipartite biadjacency matrix
    ## input
    ## Z - interaction matrix
    ## Z_est = posterior interaction matrix (optional)
    ## type  = hosts, parasites, both
    ## extra input for colours 

    ## Optional estimated presence/absence matrix (Z_est) can be added to existing plot.
    para_degrees <- as.data.frame(table(colSums(Z)))
    para_degrees$Var1 <- as.numeric(para_degrees$Var1)
    ## para_degrees = para_degrees[-which(para_degrees$Var1<2),]
    host_degrees <- as.data.frame(table(rowSums(Z)))
    host_degrees$Var1 <- as.numeric(host_degrees$Var1)
    ## host_degrees = host_degrees[-which(host_degrees$Var1<2),]

    xlim = c(1, max(para_degrees$Var1,host_degrees$Var1)*1.5)
    ylim = c(1, max(para_degrees$Freq,host_degrees$Freq)*1.5)

    if (!missing(Z_est)){
        para_est <- as.data.frame(table(colSums(Z_est)))
        para_est$Var1 <- as.numeric(para_est$Var1)
        ## para_est = para_est[-which(para_est$Var1<2),]
        host_est <- as.data.frame(table(rowSums(Z_est)))
        host_est$Var1 <- as.numeric(host_est$Var1)
        ## host_est = host_est[-which(host_est$Var1<2),]
        xlim = c(1, max(para_degrees$Var1,host_degrees$Var1,
            para_est$Var1, host_est$Var1)*1.5)
        ylim = c(1, max(para_degrees$Freq,host_degrees$Freq,
            para_est$Freq, host_est$Freq)*1.5)
    }
    gpch = c('+', '*')
    if(type=='parasites'){
        plot((para_degrees), type="p", col=parasite.col, pch=gpch[2], log="xy", xlim=xlim, ylim=ylim, ylab="Number of nodes", xlab="Degree", cex.lab = 1.5,cex.axis = 1.5)
        legend(xlim[2]*0.03, ylim[2]*1.4, c("Parasites"), col = parasite.col,
               pch = gpch[2], bty="n",pt.cex=1.5, cex=2)
        if(!missing(Z_est)){
            points((para_est), type="p", col=parasite.col, pch=16)
            legend(xlim[2]*0.03, ylim[2]*0.8, c("Est"), col = parasite.col,
               pch = 16, bty='n', pt.cex=1.5, cex=2)
        }
    }
    if(type=='hosts'){
        plot((host_degrees), type="p", col=host.col, pch=gpch[1], log="xy", xlim=xlim, ylim=ylim, ylab="Number of nodes", xlab="Degree", cex.lab = 1.5,cex.axis = 1.5)
        legend(xlim[2]*0.03, ylim[2]*1.4, c("Hosts"), col = host.col,
               pch = gpch[1], bty="n", pt.cex=1.5, cex=2)
        if(!missing(Z_est)){
            points((host_est), type="p", col=host.col, pch=16)
            legend(xlim[2]*0.03, ylim[2]*0.8, c("Est"), col = host.col,
                   pch = 16, bty='n',pt.cex=1.5, cex=2)
        }
    }
    if(type=='both'){
        plot((para_degrees), type="p", col=parasite.col, pch=gpch[2], log="xy", xlim=xlim, ylim=ylim, ylab="Number of nodes", xlab="Degree", cex.lab = 1.5,cex.axis = 1.5)
        points((host_degrees), type="p", col=host.col, pch=gpch[1])
    legend(xlim[2]*0.03, ylim[2]*1.4, c("Parasites", "Hosts"), col = c(parasite.col, host.col),
           pch = gpch[2:1], bty = 'n', pt.cex=1.5, cex=2)
        if (!missing(Z_est)) {
            points((para_est), type="p", col=parasite.col, pch=16)
            points((host_est), type="p", col=host.col, pch=16)
            legend(xlim[2]*0.03, ylim[2]*0.33, c("Est"), col = c("black"),
                   pch = 16, bty='n', pt.cex=1.5, cex=2)
        }
    }
}

sample_parameter<-function(param, MODEL,Z, tree, size = 1000, weights=NULL){
    sample_mcmc<-function(mcmc_sample,nObs,  size =1000){
        if(is.matrix(mcmc_sample))
            mcmc_sample[, sample.int(nObs, size, replace = TRUE)] else
        mcmc_sample[sample.int(nObs, size, replace = TRUE)]
    }
    if(!is.null(weights) & length(weights)!=size)
        stop('weights are not the same sampling size.')
    t.max = get.max.depth(tree)
    dist.original = cophenetic(tree)
      if(grepl('dist', MODEL)) {
        Y = W = 1
    }else{
        Y = sample_mcmc(param$y, ncol(param$y), size)
        W = sample_mcmc(param$w, ncol(param$w), size)
    }
    if(grepl('(dist|full)', MODEL)){
        Eta = sample_mcmc(param$eta, length(param$eta), size)
    }else Eta = 1

    
    zeroZ = which(Z>0)
    P <- 0
    for(s in 1:size){
        ## aux = sapply(1:size, function(s){
        ## setting affinity to 1 in distance model
        if(grepl('(aff|full)', MODEL)){
            YW = outer(Y[,s], W[,s])
        } else YW = 1
        ## Full or distance model
        ## Creating distance
        if(grepl('(full|dist)', MODEL)){
            distance = 1/EB.distance(dist.original, t.max, Eta[s])
            diag(distance)<-0
            distance = distance %*% Z
            distance[distance==0] <- if(grepl('dist', MODEL)) Inf else 1
        } else distance = 1
        ## models
        ## P = 1-exp(-outer(Y, W))                 # affinity model
        ## P = 1-exp(-distance)                    # distance model
        ## P = 1-exp(-YW*distance)                 # full model
        ## All in one probability matrix
        if(!is.null(weights)){
            a =  1-  exp(-YW*distance)
            Pg = a * weights[s] /(1-a + weights[s] * a)
            Pg[zeroZ] <- a[zeroZ]
        }else{
            P = P + 1-  exp(-YW*distance)
        }
    }
    ## })
    matrix(P/size, nrow = nrow(com), ncol = ncol(com))

}
