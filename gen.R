#!/bin/R
##################
tol.err=1e-5
## Gibbs sampling based on sections 2.4 and 2.5
rExp<-function(l,a=1){
    ## Sampling from a truncated exponential distribution
    ##l = outer(c(1,2),c(2,3))
	unif = runif(length(l))
	-log(1- unif*(1-exp(-l*a)))/(l + tol.err)
}

rExp.mean<-function(l){
    ## mean of a one-truncated exponential distribution
    (1 - (l+1)*exp(-l))/(l*(1-exp(-l)))
}

raffinity.MH<-function(old, z, Ud, sig=0.1, hyper){
    a = hyper[1]
    b = hyper[2]
    epsilon = rnorm(length(old), 0, sd = sig)
    new = old +  -2*(old + epsilon <=0)*epsilon + epsilon    

    likeli = (z+a-1)*log(new/old) -  (new - old)*(b +Ud)
    
    ratio = pmin(1, exp(likeli))
    u = 1*(runif(length(old))<=ratio)
    u*new + (1-u)*old
}

rEta.cophenetic<-function(eta.old,tree,tree.ht,pdist.old, pd0,i, Z, ywU, eta_sd =0.01){
    for(l in 1:1){
        eta.prop = eta.old + rnorm(length(eta.old), 0, sd = eta_sd)
        ##eta.prop = sign(eta.prop)*min(abs(eta.prop),20)
        dist = 1/cophenetic(eb.phylo(tree, tree.ht, eta.prop))
        diag(dist)<-0
        pdist.new = c(dist[i,]%*%Z)
        no0 = which(pd0)
        if(length(no0)){
            if(length(n0)==sum(Z[i,])) likeli = -Inf else 
            likeli = sum(((log(pdist.new)- log(pdist.old))*Z[i,])[-no0])-
                sum((ywU*(pdist.new - pdist.old))[-no0]) 
        }else{
            likeli = sum((log(pdist.new)- log(pdist.old))*Z[i,] )-
                sum(ywU*(pdist.new - pdist.old))
        }
        ratio = min(1, exp(likeli));ratio
        u = (runif(1)<=ratio)
        if(u) { eta.old  = eta.prop; pdist.old = pdist.new}
    }
    list (eta=eta.old, dist=pdist.old)
}


rEta.Horiz<-function(eta.old,dist,pdist.old,i, Z, y, w, U, eta_sd =0.01, alpha){
    ## A function that generates the prior for the power of distance
    for(l in 1:1){
        eta.prop = eta.old*exp(rnorm(1, 0, sd = eta_sd))
        eta.prop =min(eta.prop,50)
        pdist.new = exp(dist[i,]*eta.prop) %*% Z
        pdist.new[pdist.new==0]<-alpha
        
        likeli = sum((log(pdist.new)-log(pdist.old))*Z[i,] )-
            sum(U*y*w*(pdist.new - pdist.old))
        
        ratio = min(1, exp(likeli));ratio
        u = (runif(1)<=ratio)
        if(u) { eta.old  = eta.prop; pdist.old = c(pdist.new)}
    }
    list (eta=eta.old, dist=pdist.old)
}


rHyper.Horiz<-function(old,Psi, m, sig = 0.1){
    new = old*exp(rnorm(length(old), 0, sd = sig ))
    new[2]<-1
    r = sum(new[1]*log(new[2]) - lgamma(new[1]) + 
                lgamma(new[1]+m) - (new[1]+m)*log(new[2]+Psi)) -
                    sum(old[1]*log(old[2]) - lgamma(old[1]) + 
                            lgamma(old[1]+ m ) - (old[1]+m)*log(old[2] + Psi))
    ratio = min(1, exp(r));ratio
    u = 1*(runif(1)<=ratio)
    u*new + (1-u)*old
}

ICM_est<-function(Z, slice = 10, eta, distOnly=FALSE, ...){

    ## Warnings
    if(missing(Z)) stop('Interaction matrix is missing!')

    ## Setting up parameters
	n_w = ncol(Z);n_y = nrow(Z)
    n = n_w*n_y

    ## Hyper parameters
    el <-list(...)
    throw.out = if(is.null(el$burn)) floor(0.5*slice) else floor(el$burn*slice)
    a_w = if(!is.null(el$a_w)) el$a_w else 1
    a_y = if(!is.null(el$a_w)) el$a_y else 1
    b_w = b_y = 1;
    tree = if(!is.null(el$tree)) el$tree else stop('no tree passed')
    y = if(!is.null(el$y)) el$y else 1
    w = if(!is.null(el$w)) el$w else 1
    
    y_sd = if(!is.null(el$y_sd)) el$y_sd else 1
    w_sd = if(!is.null(el$wd_sd)) el$wd_sd else apply(Z, 2,sd)
    
    eta_sd = if(!is.null(el$eta_sd)) el$eta_sd else 0.1
    beta = if(!is.null(el$beta)) el$beta else 1

    print(sprintf("Run for %i slices, and %i burn ins",slice, throw.out))
    print(sprintf("Using affinity parameters: %s",!distOnly))
    print(sprintf("Using Similarity matrix: %s",!is.null(el$tree)))
    print(sprintf("Updating eta: %s",!missing(eta)))
    print(sprintf('Input matrix is Binary: %s, with dimension:', all(range(Z)==c(0,1))))
    print(dim(Z))

    ## Variable holders
    y0<-matrix(y, nrow= n_y, ncol =slice+1)
    w0<-matrix(w, nrow= n_w, ncol =slice+1)

    subItra = n_y 
    w0in<-matrix(w, nrow= n_w, ncol =subItra+1)
    y0in<-matrix(y, nrow= n_y, ncol =subItra+1)
    petain = rep(0, subItra+1)
    
	Z0 <- Z==0
    U0 <- Z*0
	mc <- colSums(1*(Z>0))
	mr <- rowSums(1*(Z>0))
    peta = rep(0, slice)

    if(is.null(el$tree)){
        pdist<- matrix(1, nrow(Z), ncol(Z))
    }
    else{
        tree.ht = arrange.tree(tree)
        dist = cophenetic(eb.phylo(tree, tree.ht, peta[1]))
        dist = 1/dist
        diag(dist)<-0
        pdist = dist%*%Z
    }
    pdist0 = pdist==0
    pdist00 = which(pdist0, arr.ind=TRUE)
    n0= NROW(pdist00)
    i=1
    s=1
    tryCatch(
        for(s in 1:slice){
            if(s%%1000==0)
                print(sprintf('Slice %d', s))
            petain[1] = peta[s]
            if(!is.null(el$tree)){
                dist = cophenetic(eb.phylo(tree, tree.ht, peta[s]))
                dist = 1/dist
                diag(dist)<-0
                pdist = dist%*%Z
                pdist[pdist00]<-1
            }
            
            ## Updatting latent scores
            yw = outer(y0[,s],w0[, s])
            U0 <- rExp(pdist*yw)
            U0[Z0] <- 1
            
            Upd = U0*pdist
            for (i in 1:subItra){
                if(!distOnly){                
                    ## Updating the parasite parameters
                    w0in[, i+1]<-raffinity.MH(w0in[,i],Z[i,],
                                              y0in[i,i]*(Upd[i,]),
                                              sig=w_sd, c(a_w, 1))
                    ## Updating host parameters
                    y0in[, i+1]<-raffinity.MH(y0in[,i],mr,
                                              (Upd)%*%w0in[,i+1],
                                              sig=y_sd, c(a_y, 1))
                }
                ## Updating similarity matix parameter
                if(!missing(eta) & !is.null(el$tree)){
                    new.eta = rEta.cophenetic(petain[i],tree,tree.ht,
                        pdist[i,],pdist0[i,],i,Z, (w0[,s]*y0in[i,i+1]*U0[i,]),
                        eta_sd)
                    petain[i+1] = new.eta$eta
                    pdist[i,] = new.eta$dist
                }
            }
            
            ## MH Adaptiveness
            # Parasite parameters (w)
            ac =1- rowMeans(abs(w0in[,1:subItra] - w0in[,1:subItra+1])<tol.err)
            w_sd = w_sd*exp(beta*(ac-0.44)/log(1 +s))
            
            ac =1- rowMeans(abs(y0in[,1:subItra] - y0in[,1:subItra+1])<tol.err)
            y_sd = y_sd*exp(beta*(ac-0.44)/log(1 +s))

            ## Tree scaling parameter (eta)
            ac =1- mean(abs(petain[1:subItra] - petain[1:subItra +1])<tol.err)
            eta_sd = eta_sd*exp(beta*(ac-0.44)/log(1 +s))
            
            if(!distOnly){
                ## w0[,s+1]<- rowSums(w0in[,-1])/subItra
                w0[,s+1] <- colSums(t(w0in[,-1])*mr)/sum(mr)
                w0in<-matrix(w0[,s+1], nrow= n_w, ncol =subItra+1)

                y0[,s+1]<- rowSums(y0in[,-1])/subItra
                y0in<-matrix(y0[,s+1], nrow= n_y, ncol =subItra+1)
                
            }
            peta[s+1]<- sum(petain[-1]*mr)/sum(mr)
            petain = rep(peta[s+1], subItra+1)
        }
       ,warning = function(w)
           {print(c('warning at (s,i):', c(s,i))); print(w);traceback()},
        error =function(e)
            {print(c('error at (s,i):', c(s,i))); print(e);traceback()} ,
        finally = print("Done!"))
    ## throwing out the burn_in stage.(approx 30%)
    if(throw.out==0) throw.out = c(1:slice) else throw.out = -c(1:throw.out)
    y0 =  y0[,throw.out] 
    w0 =  w0[,throw.out]
    if(!missing(eta)){
        eta = peta[throw.out]
    }else eta = NULL
    param = list(w = w0, y = y0, burn_in = slice - max(-throw.out), throw.out = max(-throw.out),eta = eta, sd = list(w=w_sd, y = y_sd, eta= eta_sd))
    param
}
##################################################
##################################################
## Tree function

arrange.tree<-function(phy){
    ## taken from .b.phylo from
    ## https://github.com/mwpennell/geiger-v2/blob/master/R/utilities-phylo.R#L1353-L1382
    ht=heights.phylo(tree)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), phy$edge[,2])
	ht$t1=Tmax-ht$end[phy$edge[mm,1]]
	ht$t2=ht$start-ht$end+ht$t1
    ht
}

eb.phylo=function(phy, heights, a){
    ## modified form .eb.phylo
    ## https://github.com/mwpennell/geiger-v2/blob/master/R/utilities-phylo.R
	## #ht=heights.phylo(phy)
	## N=Ntip(phy)
	## Tmax=ht$start[N+1]
	## mm=match(1:nrow(ht), phy$edge[,2])
	## ht$t1=Tmax-ht$end[phy$edge[mm,1]]
	## ht$t2=ht$start-ht$end+ht$t1
    if(a==0) return(phy)
    bl = (exp(a*heights$t2)-exp(a*heights$t1))/(a)
    phy$edge.length=bl[phy$edge[,2]]
    phy
}


heights.phylo=function(x){
    phy=x
	phy <- reorder(phy, "postorder")
	n <- length(phy$tip.label)
	n.node <- phy$Nnode
	xx <- numeric(n + n.node)
	for (i in nrow(phy$edge):1) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
	root = ifelse(is.null(phy$root.edge), 0, phy$root.edge)
	labs = c(phy$tip.label, phy$node.label)
	depth = max(xx)
	tt = depth - xx
	idx = 1:length(tt)
	dd = phy$edge.length[idx]
	mm = match(1:length(tt), c(phy$edge[, 2], Ntip(phy) + 1))
	dd = c(phy$edge.length, root)[mm]
	ss = tt + dd
	res = cbind(ss, tt)
	rownames(res) = idx
	colnames(res) = c("start", "end")
	res = data.frame(res)
	res
}


