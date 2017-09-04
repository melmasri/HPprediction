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

	image(1:ncol(Z), 1:nrow(Z), t(Z[nrow(Z):1,]),
		col = c('white', 'black'), ylab=ylab, xlab=xlab,
		useRaster=TRUE,srt=45, axes=FALSE)
	axis(1, at = 5*0:(ceiling(ncol(Z)/5)))
	axis(2, at = c(1,10*1:(ceiling(nrow(Z)/10))), labels = c(10*(ceiling(nrow(Z)/10)):1,1))
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
    require(reshape)
    P[Z>0]<--1
    aux =   melt(P)
    aux = aux[order(aux$value,decreasing=TRUE),]
    colnames(aux)<-c('Hosts', 'Parasites', 'p')
    aux[1:topX,]
}


cross.validate.fold<-function(Z, n= 10, min.per.col = 1){
    ## n-fold cross validation
    ## Returns a matrix of 3 columns, the first two are the (row,col) index of the pair,
    ## the third is the group
    if(max(range(Z))>1) Z[Z>0]<-1
    pairs = which(Z==1, arr.ind=T)
    
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
            b = a[sample(length(a),1)]
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

