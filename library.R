#!/bin/R
##libraries
## Useful function
 `%+%` <- function(a, b) paste0(a, b)
 `% %` <- function(a, b) paste(a, b)

### supporting function
# A CRM Lambda(dw , dtheta) = lambda(w)*h(theta) dw dtheta
# h is the location density
# lambda is the mass density
# h is independent from lambda and is an i.i.d.

## Library of functions
lambda<-function(w){
	# lambda is the Levy intensity of the generalize gamma process GGP
	# Generalized Gamma process (GGP) (sigma =0)
	# inverse Gaussian process (IG) ( sigma =0.5)
	# stable process (SP)(tau = 0)
	# General form: (# Default is stable GGP (sigma =0, tau=0)
	(a/gamma(1-sigma))*w^(-sigma -1)*exp(-w*tau)
	}


kappa<-function(n,z, process='gamma'){
	# General form
	#kappa(n,z) = int_0^\infty lambda(w) w^n exp(-zw)dw
	# for lambda(w) = a/gamma(1-sigma) w^(-sigma -1)exp(-w tau)
	# After integration kappa(n,z) = (a/(tau +z)^(n-sigma))(gamma(n-sigma)/gamma(1-sigma))
	switch(process, gamma 		= (a/(z+tau)^n)*(gamma(n)),
					gen.gamma	= (a/(z+tau)^(n-sigma))*gamma(n-sigma)/gamma(1-sigma),
					inv.gamma	= (a/(z+tau)^(n-0.5))*gamma(n-0.5)/sqrt(pi),
					stable		= (a/z^(n-sigma))*gamma(n-sigma)/gamma(1-sigma)
		 )
}

psi<-function(t,b=0, process='gamma'){
	# In general $\psi(t) = (a*tau^\simga/gamma(1-simga)) \sum_{k=1}^\infty (-1)^(k+1) (t/tau)^k gamma(k-sigma)/k!
	# for tau =1 and sigma =0, \psi(t) reduced to a*log(1+t/tau) for |t/tau|<1. Using Taylor expansion
	# See Supplementary Material of Caron, F. (2012) Bayesian nonparametric models for bipartite graphics.off
	tau_ = tau +b
	# Howeve, the setting is augmented to include \tilde{psi}. Hence, tau becomes tau +b with default b=0
	switch(process, gamma 		= a*log(1+t/tau_ ),
					gen.gamma	= (a/sigma)*((t+tau_)^sigma - tau_^sigma)	,
					inv.gamma	= 2*a*(sqrt(t+tau_) -sqrt(tau_)),
					stable		= a*t^sigma/sigma)
}

sample_u<-function(m,m_ind, s,i){
	# Sampling from the distribution of u_(n+1)|z_(n+1,j), u_(1:n,j) as in Equation (15)
	# This is for the GGP and is a direct sampling.
	# Based on section E in the supplementary martial

	if(length(m_ind)==0) return(0)
	if(i ==1 ) stop('Sampling start by the second reader')

	# F^(-1)(y) is:
	#i =2 	# i is the current row index, for the first row a starting point is needed
	#m_ind	# column index of the active column
	#m 		# column sum of the active ones, i.e. m_j
	#s = if(i==2) y[1:(i-1)]* U[1:(i-1),m_ind]	 else y[1:(i-1)]%*% U[1:(i-1),m_ind]
	p = runif(length(m_ind))

	u_new = sapply(1:length(m), function(j){
		if(m[j]- sigma!=0){
			(1/y[i])*(((tau+s[j])^(-m[j]+sigma)+
			+ p[j]*(( tau+s[j]+y[i])^(-m[j]+sigma) -(tau +s[j])^(-m[j]+sigma)))^(1/(-m[j]+sigma)) -
			+(tau + s[j]))
		}else{
			((tau + s[j])/(y[i]))*((1+  y[i]/(tau+ s[j] ))^p[j] -1)
		}
		})
	if(any(abs(u_new)>1)) stop('Probability for u_new out of range')
	u_new
}

gen_process<-function(nodes_w, nodes_y, y){
	# See section 2.4
	U <-matrix(1, ncol=nodes_w, nrow=nodes_y)  # latent variable, should be replaced with the original scoring.
	Z<-matrix(0, ncol=nodes_w, nrow=nodes_y)  # the binary matrix
	#insuring that the first reader read at least one book.
	u_Z = runif(1)
	repeat{
		no_y = rpois(1,psi(y[1], process = 'gamma'))
		if(no_y>1) break
	}

	Z[1,1:no_y]<-1
	U[1,1:no_y]<- runif(no_y)
	K = 1:no_y 	# Bank of already counted books
	s = 	y[1] * U[1,K]  # the first s_nj, i.e., s_1j
	m = 	rep(1,no_y)
	for(i in 2:nodes_y){
		# Sampling of Z
		# probability to select from already sampled books
		prob_old = 1 - kappa(m, tau + y[i] + s, process)/kappa(m, tau + s, process)
		Z[i,K]	 <- 1*(u_Z<=prob_old)
		#Z[i,K]	 <- 1*(runif(K)<=prob_old)
		new_books = rpois(1,psi(t=y[i],b = sum(y[1:(i-1)]) , process))
		if(new_books!=0){
			Z[i,max(K) + 1:new_books] <-1
			K = 1:(max(K)+ new_books)	# Setting up also for the next i
		}

		# Sampling of U and X
		m_ind = which(Z[i,K]==1, arr.ind=TRUE)
		m = colSums(Z[1:i,K])		# This m is also set up for the next i
		s_aux = if(i==2) y[1:(i-1)]* U[1:(i-1),m_ind]	 else y[1:(i-1)]%*% U[1:(i-1),m_ind]
		new_u = sample_u( m=m[m_ind],m_ind,s_aux, i )
		U[i,m_ind]<-new_u

		## Setting up for the new i
		s = y[1:i] %*% U[1:i,K]
	}

	#Z <- Z[,K]
	#U <- U[,K]
	list(U = U, Z=Z, X = -log(U))
}


plot_Z<-function(Z, xlab, ylab, ...){
	# UseRaster for rendering large matrices more efficiently
	# If rows are readers and columns are books, the Z matrix has to be flipped and transposed
	# a t(Z[nrow(Z):1,])
	# Axis to draw the proper limits and to rotated numbering for the y-axis.
	if(missing(ylab)) ylab = 'readers'
	if(missing(xlab)) xlab = 'books'

	image(1:ncol(Z), 1:nrow(Z), t(Z[nrow(Z):1,]),
		col = c('white', 'black'), ylab=ylab, xlab=xlab,
		useRaster=TRUE,srt=45, axes=FALSE)
	axis(1, at = 5*0:(ceiling(ncol(Z)/5)))
	axis(2, at = c(1,10*1:(ceiling(nrow(Z)/10))), labels = c(10*(ceiling(nrow(Z)/10)):1,1))
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


lof_freq<-function(X, col=TRUE){
    ## Left ordering based on column rank or row rank.
    if(col)
        X[,order(colSums(X), decreasing=T)]
    else    X[order(rowSums(X), decreasing=F),]
}

lof<-function(Z, indeces = FALSE){
    ## Given a binary matrix Z. Where the rows is fixed
    ## A function that left orders the matrix sequentially from row 1 to n
    ## based on firstappearance of columns.
    if(min(range(Z))<0) stop('Range is less that 0.')
   	## if(max(range(Z))>1){
    ##     Z1 =Z
    ##     Z[Z>0]<-1
    ## }
    active_col <- apply(Z,1,function(r) which(as.vector(r)>0))
	bank = active_col[[1]]
	for(i in 1:nrow(Z)){
			a = setdiff(active_col[[i]], bank)
			bank=c(bank,a)
		}
    if(indeces) bank else  Z[,bank]
}


threshold_point <-function(param, Z, dist){
	## A function that returns The
	## thres 	: the point of intersection of the ROC curves
	## P 		: the probability matrix (graphon)
	## no.1r1 	: the ROC green curve
	## no.0r0 	: the ROC red curve
	p_burn = param$burn_in +1 - 1:floor(0.3*param$burn_in)
	hw = rowMeans(param$w[, p_burn])
	hy = rowMeans(param$y[, p_burn])
	if(missing(dist)){
		hP = 1-exp(-outer(hy, hw))
	}else {
		hP= 1-exp(-(outer(hy,hw)*dist))
	}
	u_seq <-seq(0,1, by = 0.01)
	res = sapply(u_seq, function(u){
			hZ = (u<=hP)*1
			c(sum((abs(hZ - Z))[Z==1]),sum((abs(hZ - Z))[Z==0]))
		})
	res= t(res)

	no.ones_are_ones = 1-  res[,1]/sum(Z)
    no.zeros_are_ones = 1- res[,2]/sum(Z==0)
	thres = u_seq[which.min(abs(no.ones_are_ones - no.zeros_are_ones))]
	list(thres = thres, P = hP, no.1r1 = no.ones_are_ones,
		no.0r0= no.zeros_are_ones)
}


dist_ordering <- function(dist, com){
    ## Validating names of hosts in dist and com
    ## Testing all names in dist exist in com
    if(! all(sapply(rownames(dist), function(r) r %in% rownames(com))))
		stop('Not all hosts in dist exist in com.')
    ## Testing all names in com exist in dist
    if(! all(sapply(rownames(com), function(r) r %in% rownames(dist))))
		stop('Not all com in dist exist in dist.')
    ## Testing that names exist only once
    aux <- sapply(rownames(com), function(r) which(r==rownames(dist)))
    if(!all(sapply(aux, length)==1))
        stop('Some names in dist do not exist only once in com.')
    ## Ordering dist
    dist = dist[aux,]
    dist = dist[,aux]
    dist
}

new.dist<-function(dist, lambda.old, lambda.new){
    ## A function that returns a new distance matrix based on new lambda
    ## sim is aht similarity matrix, is ia assumed that the diag is 1.
    sim = 1/dist
    diag(sim)<-0
    d.max  = max(sim)                   # Maximum distance
    sim.new = d.max + (sim - d.max)*lambda.new/lambda.old
    diag(sim.new)<-1
    1/sim.new
}

get.P.mode <-function(param,Z, dist, lambda.old){
    ## A function that returns the P matrix with the mode of the variables
    mod.var = getMode(param)
    P = outer(mod.var$y, mod.var$w)
    ## Integrating the dist
    if(!missing(dist) && !missing(Z)){
        ## If beta is available
        if(any(grepl('beta', names(mod.var), perl=T)))
            P = P*mod.var$beta

        if(any(grepl('lambda', names(mod.var), perl=T)) && !missing(lambda.old)){
            P = P*(new.dist(dist, lambda.old, mod.var$dist_lambda)%*%Z)
        }else{
            P = P*(dist%*%Z)
            if(any(grepl('lambda', names(mod.var), perl=T)) && missing(lambda.old))
                warning('Mode of lambda is available but not used because lambda.old is missing')

        }
    }

    1 - exp(-P)
}

rocCurves<-function(param,Z,Z_cross, dist, lambda.old,plot=TRUE,P,bins=200, all=FALSE){
    ## Computes the ROC curves as x = FPR and y = TPR
    ##  FPR = FP/(FP +TN)
    ##  TPR = TP/(TP+FN)
    if(missing(P)){
        if(missing(Z_cross))
            P = get.P.mode(param, Z, dist, lambda.old) else
    P = get.P.mode(param, Z_cross, dist,lambda.old)
    }

    if(length(unique(as.vector(Z)))!=2) Z = 1*(Z>0)
    if(length(unique(as.vector(Z_cross)))!=2) Z_cross = 1*(Z_cross>0)
    u = seq(0, 1, length.out = bins)
    Z1 = (Z - Z_cross)==1
    if(all) Z1 = Z==1
    m = sum(1*(Z1))
    Z2 = Z==0
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
    roc = data.frame(u=u, t(aux), FPR = FPR, TPR = TPR)
    max.point =   which.min(abs(roc$TPR-1) + roc$FPR)
    threshold = roc$u[max.point]
    n = length(u)
    auc =  0.5*t(abs(FPR[2:n] - FPR[2:n -1]))%*% (TPR[2:n] + TPR[2:n -1])
    auc = round(100*auc,2)
    if(plot){
        plot(FPR, TPR, xlab='FPR (1-specifity)', ylab = 'TPR (sensitivity)',
             type ='l', col='red', lwd=2, main = paste('AUC ', auc), xlim = c(0,1), ylim = c(0,1), pch =6, lty=4)
        abline(a = 0, b=1,col='black',lty=2, lwd=2)
    }

    list(roc= roc,auc = auc, threshold = threshold,
         max.point= c(FPR[max.point], TPR[max.point]) , P= P)
}


getMode<-function(param){
    ## A function that returns the posterior mode of all paramters of the model
    extra.throw.out = 1
    if(param$throw.out/(param$throw.out + param$burn_in) <=0.1){
        extra.throw.out =1:round(0.2*(param$burn_in + param$throw.out) - param$throw.out)
    }
    estimate_mode <- function(x) {
        d <- density(x[-extra.throw.out])
        d$x[which.max(d$y)]
    }
    names =  names(param)[!sapply(param, is.null)]
    names = names(param)[-grep('(etashoot|throw.out|burn_in|w_star|Z)',names(param))]
    
    aux = lapply(param[names], function(r){
        if (!is.null(dim(r)))
                              apply(r,1,estimate_mode) else
                          if(length(r)>1) estimate_mode(r)
                      })
    aux
}

getMean<-function(param){
    ## A function that returns the posterior mode of all paramters of the model
    extra.throw.out = 1
    if(param$throw.out/(param$throw.out + param$burn_in) <=0.1){
        extra.throw.out =1:round(0.2*(param$burn_in + param$throw.out) - param$throw.out)
    }
    estimate_mode <- function(x)    mean(x)

    names =  names(param)[!sapply(param, is.null)]
    names = names(param)[-grep('(hh|sd|etashoot|throw.out|burn_in|w_star|Z)',names(param))]
    
    aux = lapply(param[names], function(r){
        if (!is.null(dim(r)))
            apply(r,1,estimate_mode) else
        if(length(r)>1) estimate_mode(r)
    })
    aux
}

topPairs<-function(P,Z,topX=20){
   ## Plotting the top new pairs
    ## Note that P is a probability matrix, not a binary one.
    require(reshape2)
    Z = com_pa
    P[Z==1]<--1
    aux=   melt(P)
    aux = aux[order(aux$value,decreasing=TRUE),]
    colnames(aux)<-c('Hosts', 'Parasites', 'p')
    aux[1:topX,]
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

cross.validate.fold<-function(Z, n= 10){
    ##  n-fold cross validation
    ## Returns a matrix of 3 columns, the first two are the (row,col) index of the pair,
    ## the third is the group
    if(max(range(Z))>1) Z[Z>0]<-1
    pairs = which(Z==1, arr.ind=T)
    col1 = which(colSums(Z)<3)
    col2 = which(colSums(Z)==3)
    if(length(col1)>0)
        pairs = pairs[-which(pairs[,'col'] %in% col1),]
    size = floor(nrow(pairs)/n)
    gr = rep(1:n, each = size)
    if(nrow(pairs) %% size!=0)
        gr = c(gr, rep(n,(nrow(pairs) %% size) ))
    pairs = cbind(pairs[sample(nrow(pairs), nrow(pairs)), ], gr = gr)

    col.good = setdiff(unique(pairs[,'col']), col2)
    
    aux = lapply(as.vector(col2), function(r){
        s =as.vector(which(pairs[,'col'] == r))
        same  = pairs[s, 'gr']
        l = same
        if(same[2]==same[1])
            l =  c(same[1],sample((1:n)[-same[1]],1), same[3])
        if(same[3]==same[1])
            l =  c(l[1],l[2], sample((1:n)[-l[1:2]],1))
        if(l[3]==l[2])
            l =  c(l[1],l[2], sample((1:n)[-l[1:2]],1))
        cbind(s, l) 
    })
    
    aux = do.call('rbind', aux)
    pairs[aux[,1], 'gr']<-aux[,2]
    print(sprintf("Actual cross-validation rate is %0.3f" , table(pairs[,'gr'])/sum(1*(Z>0))))
    pairs[order(pairs[,'gr']),]
    
}



postertior.predictive<-function(Z,Z_cross, param, dist, lambda.old){
    ## A function that calculates the posterior predictive
    ## given a list of the w and y parameters of the gibbs sampling
    ##  in the log scale, for a single simulation of n_1 gibbs samples
    ## the posterior is
    ## (1/n_1)* sum(log(f_k(Z|w^j, y^j, D_k^j)) )
    ## where f_k is the likelihood given the other parameters
    n = param$burn_in
    tol.err = 1e-16
    mod.var = getMode(param)
    Z1 = Z_cross==1
    if(!missing(dist)){
        if(grepl('(?=.*beta)(?=.*lambda)', paste0(names(param), collapse='-'), perl=TRUE)){
            ## Both beta and lambda are available
            if(missing(lambda.old))
                stop('Gibbs for lambda_phy is available but no original lambda')
            post =  sapply(1:n, function(r){
                               newdist = new.dist(dist, lambda.old, param$dist_lambda[r])%*%Z_cross
                               P = 1- exp(-outer(param$y[,r], param$w[,r])*param$beta[r]*newdist)
                               P[Z1]<-1 - tol.err
                               P[P==0]<-tol.err
                               sum(Z*log(P) + (1-Z)*log(1-P))
                           })

        }else if(any(grepl('beta', names(param)))){
             ## Only beta is available
            newdist = dist%*%Z_cross
              post =  sapply(1:n, function(r){
                  P = 1-exp(-outer(param$y[,r], param$w[,r])*newdist*param$beta[r])
                  P[Z1]<-1 - tol.err
                  P[P==0]<-tol.err
                  sum(Z*log(P) + (1-Z)*log(1-P))
              })
        }else if(any(grepl('lambda', names(param)))){
            ## only Lambda is available
            post =  sapply(1:n, function(r){
                newdist = new.dist(dist, lambda.old, param$dist_lambda[r])%*%Z_cross
                P = 1- exp(-outer(param$y[,r], param$w[,r])*newdist)
                P[Z1]<-1 - tol.err
                P[P==0]<-tol.err
                sum(Z*log(P) + (1-Z)*log(1-P))
            })

        }else{
            ## Only dist is available
            newdist = dist%*%Z_cross
            post =  sapply(1:n, function(r){
                P = 1-exp(-outer(param$y[,r], param$w[,r])*newdist)
                P[Z1]<-1 - tol.err
                P[P==0]<-tol.err
                sum(Z*log(P) + (1-Z)*log(1-P))
            })
        }
    }else{
        ## Dist is not available only w, and y
        post =  sapply(1:n, function(r){
            P = 1- exp(-outer(param$y[,r], param$w[,r]))
            P[Z_cross==1]<-1 - tol.err
            P[P==0]<-tol.err
            sum(Z*log(P) + (1-Z)*log(1-P))
        })
    }
    mean(post)
}

get.P.all<-function(Z,Z_cross, param, dist, lambda.old){
    ## A function that retuns a list of all probability matrices
    n = param$burn_in
    tol.err = 1e-16
    if(!missing(dist)){
        if(grepl('(?=.*beta)(?=.*lambda)', paste0(names(param), collapse='-'), perl=TRUE)){
            ## Both beta and lambda are available
            if(missing(lambda.old))
                stop('Gibbs for lambda_phy is available but no original lambda')
            post =  lapply(1:n, function(r){
                               newdist = new.dist(dist, lambda.old, param$dist_lambda[r])%*%Z_cross
                               1- exp(-outer(param$y[,r], param$w[,r])*param$beta[r]*newdist)
                           })
        }else if(any(grepl('beta', names(param)))){
             ## Only beta is available
            newdist = dist%*%Z_cross
            post =  lapply(1:n, function(r)
                1-exp(-outer(param$y[,r], param$w[,r])*newdist*param$beta[r]) )
        }else if(any(grepl('lambda', names(param)))){
            ## only Lambda is available
            post =  lapply(1:n, function(r){
                newdist = new.dist(dist, lambda.old, param$dist_lambda[r])%*%Z_cross
                1- exp(-outer(param$y[,r], param$w[,r])*newdist)
            })
        }else{
            ## Only dist is available
             newdist = dist%*%Z_cross
             post =  lapply(1:n, function(r) 1-exp(-outer(param$y[,r], param$w[,r])*newdist))
         }
    }else{
        ## Dist is not available only w, and y
         post =  lapply(1:n, function(r)  1- exp(-outer(param$y[,r], param$w[,r])))
     }
    post
}

AUC<-function(Z,Z_cross, param, dist, lambda.old){
    ## A function that returns the AUC based on the mathematical formula.

    n = min(500,param$burn_in)
    tol.err = 1e-16
    pairs  = which(Z - Z_cross ==1, arr.ind =T)
    m1= nrow(pairs)
    m = sum(1-Z)

    if(!missing(dist)){
        if(grepl('(?=.*beta)(?=.*lambda)', paste0(names(param), collapse='-'), perl=TRUE)){
            ## Both beta and lambda are available
            if(missing(lambda.old))
                stop('Gibbs for lambda_phy is available but no original lambda')
            post = lapply(1:n, function(r){
                              newdist = new.dist(dist, lambda.old, param$dist_lambda[r])%*%Z_cross
                              P = 1- exp(-outer(param$y[,r], param$w[,r])*param$beta[r]*newdist)
                              sum(apply(pairs, 1, function(x)  sum((P[x[1],x[2]]> P)*1*(1-Z))))
                          })
        }else if(any(grepl('beta', names(param)))){
             ## Only beta is available
             newdist = dist%*%Z_cross
             post =  lapply(1:n, function(r){
                                P = 1-exp(-outer(param$y[,r], param$w[,r])*newdist*param$beta[r])
                                sum(apply(pairs, 1, function(x)  sum((P[x[1],x[2]]> P)*1*(1-Z))))
                            })
        }else if(any(grepl('lambda', names(param)))){
            ## only Lambda is available
            post =  lapply(1:n, function(r){
                newdist = new.dist(dist, lambda.old, param$dist_lambda[r])%*%Z_cross
                P = 1- exp(-outer(param$y[,r], param$w[,r])*newdist)
                sum(apply(pairs, 1, function(x)  sum((P[x[1],x[2]]> P)*1*(1-Z))))
            })
        }else{
            ## Only dist is available
             newdist = dist%*%Z_cross
             post =  lapply(1:n, function(r){
                                P=1-exp(-outer(param$y[,r], param$w[,r])*newdist)
                                sum(apply(pairs, 1, function(x)  sum((P[x[1],x[2]]> P)*1*(1-Z))))
                            })
         }
    }else{
         ## Dist is not available only w, and y
         post =  lapply(1:n, function(r){
                            P = 1- exp(-outer(param$y[,r], param$w[,r]))
                            sum(apply(pairs, 1, function(x)  sum((P[x[1],x[2]]> P)*1*(1-Z))))
                        })
     }
    sum(unlist(post))/(n*m1*(m-m1))
}

#### Document these function

plot_thresh<-function(i,y_est,w_est,Z,host = TRUE, ...){
	m<-match.call()
	miss_y <-is.null(m$y_o)

	if(!miss_y) y<-m$y_o
	miss_w <-is.null(m$w_o)
	if(!miss_w) w<-m$w_o
	EulerMas = 0.5772
	# Plots the threshold of of the estimation
	if(!miss_y)
		if(length(y) != length(y_est)) stop('Estimate and original parameter are not the same length.')
	if(!miss_w)
		if(length(w) != length(w_est)) stop('Estimate and original parameter are not the same length.')

	if(host){
		if(!miss_y)	thresh_o = -log(y[i])  # Original threshold
		thresh_e = -log(y_est[i])
		if(!miss_w) S_o = log(w) + EulerMas
		S_e = log(w_est) + EulerMas
	}else{
		if(!miss_w) thresh_o = -log(w[i])  # Original threshold
		thresh_e = -log(w_est[i])
		if(!miss_y) S_o = log(y) + EulerMas
		S_e = log(y_est) + EulerMas
	}
	perc =1.5
	mini   =  if(!miss_y & !miss_w) min(S_e, S_o) else min(S_e)
	maxi   = 	max(S_e, -thresh_e*perc,thresh_e*perc)
	#maxi   = 	max(S_e, S_o , -thresh_e*1.5, -thresh_o*1.5,thresh_e*1.5, thresh_o*1.5)
	n 		= length(S_e)

   	dev.new()
    # plotting estimated points
	plot(S_e, n:1, type ='p', pch = '*', xlim = c(mini, maxi),
		col='red',
		xlab = 'Scores: red-> estimated, blue -> original',
		ylab = if(!host) 'Hosts (rows)' else 'Parasite (columns)',
         main = if(host) paste('Host', i)  else paste('Parasite', i))
        # draw lines from minimum to the estimated points ( horizontal lines)
	for(j in 1:n)
		lines(c(mini, S_e[j]), c(n-j+1, n-j+1), type ='l', col='red')
    # draw estimated threshold
    abline(v =thresh_e, col ='red')
    # if original y and w are avialable plot them in blue.
    if(!miss_y & !miss_w){
        lines(S_o, n:1, type ='p', pch = '*', col='blue')
        abline(v =thresh_o, col ='blue')
    }
    # if the original bipartite matrix is available add circles.
    if(!missing(Z)){
		point = if(host) which(Z[i, ]==1) else which(Z[,i]==1)
		if(sum(point)>0)
			lines(data.frame(x = S_e[point],y = n-point+1) , type = 'p', col='black')
	}
}

generate_interactions<-function(r, c, eta = 0.3){
    ## Generate a sample of Z, D, w and y
    i = floor(r*1.5)
    j = floor(c*1.5)
    distance_matrix<-function(n){
        d = matrix(0, ncol=n, nrow=n)
        a = rexp(sum(upper.tri(d)),5)
        d[upper.tri(d)]<-a
        d = t(d)
        d[upper.tri(d)]<-a
        d
    }
    D = distance_matrix(i)
    De = D^eta
    w = rgamma(j ,0.2, 1)
    y = rgamma(i, 0.2, 1)
    Z = sapply(1:j, function(s){
        v= rep(0, i)
        ra = sample(1:i,i)
        v[ra[1]]<-1
        ## Setting of first interaction
        for(k in ra[2:length(ra)]){
            p = 1- exp(-y*(w[s]*(v%*%(De))))
            v[k]<-1*(runif(1)<=p[k])
        }
        v
    })
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
            print("No. columns less than r")
    }
    ## Removing row numbers
    aux = which(rowSums(Z)<=1)
    Z = Z[-aux,]
    y = y[-aux]
    D = D[-aux,]
    D = D[,-aux]
    ## Re ordering
    bank = lof(Z, indeces = TRUE)
    w = w[bank]
    Z = Z[,bank]
    list(w=w, y=y, eta= eta,Z=Z, D=D)
}

ana.table<-function(com, comCross, roc, plot=FALSE){
    com = 1*(com>0)
    comCross  = 1*(comCross>0)
    zz = 1*(roc$P>roc$threshold)
    if(plot) plot_Z(zz)
    tb= data.frame(auc = roc$auc, thresh= roc$threshold,
        tot.inter = sum(com), hold.out= sum(abs(comCross - com)[com==1]),
        pred = sum(zz[com==1 & comCross==0])/sum(abs(comCross - com)[com==1]),
        pred.all = sum(zz[com==1])/sum(com))
     ##   zeros =sum(comCross==0),
    ##     TP =sum(zz[com>0]),
    ##     FN =sum(1-zz[com>0]),
    ##     FP =sum(zz[com==0]),
    ##     TN = sum(1-zz[com==0]),
    ##     no.hosts = dim(com)[1],
    ##     no.parasite = dim(com)[2])
    ## 
    tb
}

ana.plot<-function(pg, com){
    r = which.max(rowSums(com));r
    c = which.max(colSums(com));c
    burn = -c(1:3)
    par(mfrow=c(3,1))
    plot(pg$w[r,burn], type='l', main = 'Most popular host')
    plot(pg$y[c,burn], type='l', main = 'Most popular parasite')
    plot(pg$eta[burn], type='l', main = 'Beta posterior')
    if(!is.null(pg$g)){
        dev.new()
        par(mfrow=c(1,1))
        hist(pg$g[burn], main = 'G posterior')
    }
}
