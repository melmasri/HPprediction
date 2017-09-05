#!/bin/R
### ==================================================
### Secondary functions
### ==================================================

getMean<-function(param){
    ## A function that returns the posterior mean of all paramters of the model
    names = names(param)[grep('(w|y|eta|g)',names(param))]

    lapply(param[names], function(r){
        if (!is.null(dim(r)))
            rowMeans(r) else
        if(length(r)>1) mean(r)
    })
    
}

plot_degree <- function(Z, Z_est, type='both', host.col='blue', parasite.col='red'){
	# Provide presence/absence matrix (Z) with rows as hosts and columns as parasites.
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
        plot((para_degrees), type="p", col=parasite.col, pch=gpch[2], log="xy", xlim=xlim, ylim=ylim, ylab="Number of Nodes", xlab="Degree")
        legend(xlim[2]*0.3, ylim[2]*0.6, c("Parasites"), col = parasite.col,
               pch = gpch[2], box.col="white")
        if(!missing(Z_est)){
            points((para_est), type="p", col=parasite.col, pch=16)
            legend(xlim[2]*0.3, ylim[2]*0.4, c("Estimated"), col = parasite.col,
               pch = 16, box.col="white")
        }
    }
    if(type=='hosts'){
        plot((host_degrees), type="p", col=host.col, pch=gpch[1], log="xy", xlim=xlim, ylim=ylim, ylab="Number of Nodes", xlab="Degree")
        legend(xlim[2]*0.3, ylim[2]*0.6, c("Hosts"), col = host.col,
               pch = gpch[1], box.col="white")
        if(!missing(Z_est)){
            points((host_est), type="p", col=host.col, pch=16)
            legend(xlim[2]*0.3, ylim[2]*0.4, c("Estimated"), col = host.col,
                   pch = 16, box.col="white")
        }
    }
    if(type=='both'){
        plot((para_degrees), type="p", col=parasite.col, pch=gpch[2], log="xy", xlim=xlim, ylim=ylim, ylab="Number of Nodes", xlab="Degree")
        points((host_degrees), type="p", col=host.col, pch=gpch[1])
    legend(xlim[2]*0.3, ylim[2]*0.6, c("Parasites", "Hosts"), col = c(parasite.col, host.col),
           pch = gpch[2:1], box.col="white")
        if (!missing(Z_est)) {
            points((para_est), type="p", col=parasite.col, pch=16)
            points((host_est), type="p", col=host.col, pch=16)
            legend(xlim[2]*0.3, ylim[2]*0.33, c("Estimated"), col = c("black"),
                   pch = 16, box.col="white")
        }
    }
}


ana.plot<-function(pg, wait = TRUE){
    rho_col <- "blue"
    gamma_col <- "red"
    eta_col <- "darkgreen"
    g_col<-'ivory4'
    
    r = which.max(rowMeans(pg$y));r
    c = which.max(rowMeans(pg$w));c
    burn = pg$burn_in - 10000:1
    burn = burn[burn>0]
    
    nrow=0
    nrow = nrow + ifelse(max(pg$y) != min(pg$y), 1, 0)
    nrow = nrow + ifelse(max(pg$w) != min(pg$w), 1, 0)
    nrow = nrow + ifelse(!is.null(pg$eta), 1, 0)
    nrow = nrow + ifelse(!is.null(pg$g), 1, 0)
    
    par(mfrow=c(nrow,1))

    if(max(pg$y) != min(pg$y)) 
        plot(pg$y[r,burn], type='l', main='', col=gamma_col, xlab='Iteration', ylab='Host parameter')

    if(max(pg$w) != min(pg$w)) 
        plot(pg$w[c,burn], type='l', main = '', col=rho_col,  xlab = 'Iteration', ylab='Parasite parameter')

    if(!is.null(pg$eta))
        plot(pg$eta[burn], type='l', main = '', col=eta_col, xlab = 'Iteration', ylab='Scaling parameter')

    if(!is.null(pg$g))
        plot(pg$g[burn], type='l', main = '', col=eta_col, xlab = 'Iteration', ylab='Uncertainty parameter')
    
    ## if(!is.null(pg$g)){
    ##     par(mfrow=c(1,1))
    ##     hist(pg$g[burn],xlab='Posterior estimate of g',main='', col=g_col)
    ## }

    if(wait)
        readline(prompt="Press [enter] to continue")

    par(mfrow=c(nrow,1))
    if(max(pg$y) != min(pg$y)){
        acf(pg$y[r,burn],lag.max=100, main='Host - most active')
        text(50, 0.8, paste0('Effective sample size: ',round(effectiveSize(pg$y[r,burn]))))
        round(effectiveSize(pg$y[r,burn]))
    }
    if(max(pg$w) != min(pg$w)) {
        acf(pg$w[c,burn],lag.max=100, main='Parasite - most active')
        text(50, 0.8, paste0('Effective sample size: ',round(effectiveSize(pg$w[c,burn]))))
    }
    if(!is.null(pg$eta)){
        acf(pg$eta[burn],lag.max=100, main='Scale parameter')
        text(50, 0.8, paste0('Effective sample size: ',round(effectiveSize(pg$eta[burn]))))
    }

    if(!is.null(pg$g)){
        acf(pg$g[burn],lag.max=100, main='Uncertainty parameter')
        text(50, 0.8, paste0('Effective sample size: ',round(effectiveSize(pg$g[burn]))))
    }
    
}

generate_interactions<-function(r, c, eta = 0.3, aj=0.5, ai=0.5, D){

    ## Generate a sample of Z, D, w and y
    i = floor(r*1.5)
    j = floor(c*2)
    print(sprintf('Attempting an %ix%i matrix, with (r,c)=(%i,%i) truncation', i,j,r,c))

    distance_matrix<-function(n){
        d = matrix(0, ncol=n, nrow=n)
        a = rexp(sum(upper.tri(d)),5)
        d[upper.tri(d)]<-a
        d = t(d)
        d[upper.tri(d)]<-a
        d
    }
    if(missing(D))
        D = distance_matrix(i)
    else{
        print(sprintf('Passed similarity matrix of dim %ix%i', nrow(D), ncol(D)))        
        if(i >= nrow(D)){
            i = nrow(D)
            warning(paste0('The number of hosts is lowered to: ', nrow(D)))
        }else{
            minset = sample(1:nrow(D), i)
            D = D[minset, ]
            D = D[, minset]
        }
    }

    De = D^eta
    w = rgamma(j ,aj, 1)
    y = rgamma(i, ai, 1)
    ## Setting of first interaction
    yw = outer(y,w)
    
    Z=matrix(0, i,j)
    print(sprintf('Running an MCMC for %i iterations',100*i))
    ##Sampling interactions
    for(s in 1:(100*i)){
        sam = sample(1:i, j, replace=TRUE)
        Diag = cbind(sam,1:j)
        dd= (De%*%Z)[Diag]
        dd = 1*(dd==0) + dd*1*(dd>0)
        Z[Diag]<- 1*(runif(j)<= 1-exp(-yw[Diag]*dd))
    }

    print('Removing single host parasites...')
    ## Removing empty or extra rows
    aux = which(colSums(Z)<=1)
    Z = Z[,-aux]
    w = w[-aux]
    if(ncol(Z)>c){
        aux = order(colSums(Z), decreasing=FALSE)
        aux = aux[1:(ncol(Z)-c)]
        Z = Z[,-aux]
        w = w[-aux]
    }else{
        if(ncol(Z)!=c)
            stop("No. columns less than c")
    }
    ## Removing row numbers
    aux = which(rowSums(Z)<1)
    if(length(aux)>0){
        Z = Z[-aux,]
        y = y[-aux]
        D = D[-aux,]
        D = D[,-aux]
    }
    ## Re ordering
    bank = lof(Z, indeces = TRUE)
    w = w[bank]
    Z = Z[,bank]
    list(w=w, y=y, eta= eta,Z=Z, D=D)
}

