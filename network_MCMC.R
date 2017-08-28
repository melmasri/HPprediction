#!/bin/R
## Tolerance error 
tol.err=1e-5

### ==================================================
### Main MCMC function
### ==================================================
network_est<-function(Z, slices = 10, tree = NULL, model.type = c('full', 'distance', 'affinity'), uncertainty = FALSE, ... ){
    require(geiger)
    ## Running options:
    ## 1 - Full to combined 2 and 3;
    ## 2 - Affinity-only model (affinity);
    ## 3 - distance-only model (distance).
    model.type= model.type[1]
    ## General warnings are checks
    if(missing(Z)) stop('Interaction matrix is missing!')
    if(!all(range(Z)==c(0,1))){
        warning('Z is converted to binary!', immediate. = TRUE, call.= FALSE)
        Z = 1*(Z>0)
    
    }
    if(any(rowSums(Z)==0)){
        stop('Z has empty rows, please remove!', immediate. = TRUE, call.= FALSE)
        Z = Z[which(rowSums(Z)>0),]
    } 
    if(any(rowSums(Z)==0)){
        stop('Z has empty columns, please remove!', immediate. = TRUE, call.= FALSE)
        Z = Z[,which(colSums(Z)>0)]
    } 
    if(slices==0)
        stop('no. of slices cannot be 0!')

    ## For affinity-only model 
    if(grepl('aff', model.type,ignore.case = TRUE)){
        if(!is.null(tree))
            warning('affinity model is chosen; ignoring tree!',
                    immediate. = TRUE, call. = FALSE)
        print('Running affinity model...')
        param = fullJoint_est(Z, iter = slices, uncertainty = uncertainty, ...)
        return (list(param=param, Z = Z))
    }

    ##  Full and distance model
    if(grepl('(dist|full)', model.type)){
        if(is.null(tree))
            stop('distance model is chosen, but tree is null!')

        if(!is.phylo(tree))
            stop('tree must be a phylogeny tree, see gieger!')

        ## testing that all tips exist in Z
        if(!all(tree$tip.label %in% rownames(Z))){
            warning('not all species in tree exist in Z; missing are removed from tree!',
                    immediate.= TRUE, call. = FALSE)
            tree <- drop.tip(tree,
                             tree$tip.label[!(tree$tip.label %in% rownames(Z))])
        }
        
        ## Testing all names in com exist in dist
        if(!all(rownames(Z) %in% tree$tip.label)){
            warning('not all row-species in Z exist tree; missing are removed from Z!',
                    immediate.= TRUE, call. = FALSE)
            Z <- Z[rownames(Z) %in% tree$tip.label,]
        }

        ## Testing that names in Z exist only once
        if(!all(sapply(sapply(rownames(Z),
                              function(r) which(r==tree$tip.label)), length)==1))
            stop('some row-species in Z exist more than once tree!')

        ## Making sure the order of hosts in Z and tree are the same
        aux = cophenetic(tree)
        row.order <- sapply(rownames(aux), function(r) which(r==rownames(Z)))
        print('Ordering the rows of Z to match tree, and left ordering the columns..')
        Z = lof(Z[row.order,])

        ## Running the MCMC
        print(paste0('Running ',
                     ifelse(grepl('dist', model.type),
                            'distance model...', 'full model...')))
        param  = ICM_est(unname(Z),
            tree = tree,slices,  distOnly = grepl('dist', model.type),
            uncertainty = uncertainty, ...)
        return(list(param = param,tree=tree, Z=Z))
    }
}

ICM_est<-function(Z, tree = tree, slices = 10, distOnly = FALSE, uncertainty = FALSE, ...){
    ## loading extra args
    el <-list(...)
    ## parameters set-up
    nw = ncol(Z);ny = nrow(Z)
    n = nw*ny                          

    y = ifelse(!is.null(el$y), el$y , 1)
    w = ifelse(!is.null(el$w), el$w , 1)
    a_w = ifelse(!is.null(el$a_w), el$a_w, 1)
    a_y = ifelse(!is.null(el$a_w), el$a_y,  1)
    b_w = b_y = 1;
    y_sd = ifelse(!is.null(el$y_sd), el$y_sd, 0.2)
    w_sd = ifelse(!is.null(el$wd_sd) ,el$wd_sd,0.2)
    eta = ifelse(is.null(el$eta), 0, el$eta)
    eta_sd = ifelse(!is.null(el$eta_sd), el$eta_sd, 0.005)

    ## Burn in set-up
    beta = ifelse(!is.null(el$beta), el$beta , 1)
    burn.in = ifelse(is.null(el$burn.in), floor(0.5*slices), floor(el$burn.in*slices))

    print(sprintf("Run for %i slices with %i burn-ins",slices, burn.in))
    ##    print(sprintf("Using uncertainty: %s",!is.null(el$g)))

    ## Variable holders
    ## outer loop
    y0<-matrix(y, nrow= ny, ncol =slices+1)
    w0<-matrix(w, nrow= nw, ncol =slices+1)
    g0<-rep(0, slices+1)
    peta = rep(eta, slices)

    ## inner loop
    subItra = ny 
    w0in = matrix(w, nrow= nw, ncol =subItra+1)
    y0in = matrix(y, nrow= ny, ncol =subItra+1)
    petain = rep(0, subItra+1)
    g0in = rep(0, subItra+1)

    Z0 <- Z==0
	mc <- colSums(Z)
	mr <- rowSums(Z)

    ## Arranging the tree
    tree.ht = arrange.tree(tree)
    dist = cophenetic(eb.phylo(tree, tree.ht, peta[1]))
    dist = 1/dist
    diag(dist)<-0
    pdist = dist%*%Z
    
    pdist0 = pdist==0
    pdist00 = which(pdist0, arr.ind=TRUE)
    n0= NROW(pdist00)

    i=1;s=1
    tryCatch(
        for(s in 1:slices){
            if(s%%100==0)
                print(sprintf('slice: %d', s))
            ## Arranging the tree
            dist = cophenetic(eb.phylo(tree, tree.ht, peta[s]))
            dist = 1/dist
            diag(dist)<-0
            pdist = dist%*%Z
            pdist[pdist00]<-1
            
            ## Updatting latent scores
            yw = outer(y0[,s],w0[, s])
            if(!uncertainty){
                U0 <- rExp(pdist*yw)
                U0[Z0] <- 1
            }else
                U0 <-rExp2(pdist*yw, g0[s], Z, Z0)
            
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
                new.eta = rEta.cophenetic(petain[i],tree,tree.ht,
                    pdist[i,],pdist0[i,],i,Z, (w0[,s]*y0in[i,i+1]*U0[i,]),
                    eta_sd)
                petain[i+1] = new.eta$eta
                pdist[i,] = new.eta$dist
                
                if(uncertainty){
                    g0in[i+1] = rg(Z[i,], l=U0[i,])
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
                w0in<-matrix(w0[,s+1], nrow= nw, ncol =subItra+1)

                y0[,s+1]<- rowSums(y0in[,-1])/subItra
                y0in<-matrix(y0[,s+1], nrow= ny, ncol =subItra+1)
                
            }
            peta[s+1]<- sum(petain[-1]*mr)/sum(mr)
            petain = rep(peta[s+1], subItra+1)

            if(uncertainty){
                g0[s+1] = sum(g0in[-1]*mr)/sum(mr)
                g0in = rep(g0[s+1], subItra+1)
            }
        }
       ,warning = function(war)
           {print(c('warning at (s,i):', c(s,i))); print(war);traceback()},
        error =function(e)
            {print(c('error at (s,i):', c(s,i))); print(e);traceback()} ,
        finally = print("Done!"))

    if(burn.in==0) burn.in = 1 else burn.in = 1:(burn.in+1)
    if(!distOnly){
        y0 =  y0[,-burn.in]
        w0 =  w0[,-burn.in]
    } else{
        y0=1
        w0=1
    }
    peta = peta[-burn.in]
    if(uncertainty)  g0 = g0[-burn.in] else g0 = NULL

    list(w = w0, y = y0, eta = peta, g = g0,
         burn.in = max(burn.in)-1,
         sd = list(w=w_sd, y = y_sd, eta= eta_sd))
}

fullJoint_est<-function(Z, iter = 10, uncertainty = FALSE, ...){
    ## loading extra args
    el <-list(...)

    ## parameters set-up
    nw = ncol(Z);ny = nrow(Z)
    n = nw*ny                          

    y = ifelse(!is.null(el$y), el$y , 1)
    w = ifelse(!is.null(el$w), el$w , 1)
    a_w = ifelse(!is.null(el$a_w), el$a_w, 1)
    a_y = ifelse(!is.null(el$a_w), el$a_y,  1)
    b_w = b_y = 1;
    y_sd = ifelse(!is.null(el$y_sd), el$y_sd, 0.2)
    w_sd = ifelse(!is.null(el$wd_sd) ,el$wd_sd,0.2)
    
    ## Burn in set-up
    batch.size = ifelse(!is.null(el$batch.size), el$batch.size, 50)
    beta = ifelse(!is.null(el$beta), el$beta , 1)
    burn.in = ifelse(is.null(el$burn.in), floor(0.5*iter), floor(el$burn.in*iter))

    print(sprintf("Run for %i iterations with %i burn-ins",iter, burn.in))
    ##    print(sprintf("Using uncertainty: %s",!is.null(el$g)))

    ## Variable holders
    y0<-matrix(y, nrow= ny, ncol =iter+1)
    w0<-matrix(w, nrow= nw, ncol =iter+1)
    g0<-rep(0, iter+1)
    peta = rep(0, iter+1)
    
    Z0 <- Z==0
	mc <- colSums(Z)
	mr <- rowSums(Z)
    s=1
    tryCatch(
        for(s in 1:iter){
            if(s%%100==0)
                print(sprintf('iteration %d', s))
            
            ## Updatting latent scores
            if(!uncertainty){
                U0 <- rExp(outer(y0[,s],w0[, s]))
                U0[Z0] <- 1
            }else
                U0 <-rExp2(outer(y0[,s],w0[, s]), g0[s], Z, Z0)
            
            ## Updating the parasite parameters
            w0[, s+1]<-raffinity.MH(w0[,s],mc,y0[,s]%*%U0,
                                    sig=w_sd, c(a_w, 1))
            ## Updating host parameters
            y0[, s+1]<-raffinity.MH(y0[,s],mr, U0%*%w0[,s+1],
                                    sig=y_sd, c(a_y, 1))
            ## Uncertain parameter sampling
            if(uncertainty)
                g0[s+1] = rg(Z, l=U0)
            
            ## MH Adaptiveness
            if(s%%batch.size==0){
                ss = s+1 - 1:batch.size
                ac =1- rowMeans(abs(w0[,ss] - w0[,ss+1])<tol.err)
                w_sd = w_sd*exp(beta*(ac-0.44)/log(1 +s/batch.size))
                
                ac =1- rowMeans(abs(y0[,ss] - y0[,ss+1])<tol.err)
                y_sd = y_sd*exp(beta*(ac-0.44)/log(1 +s/batch.size))
                
            }
        }
       ,warning = function(war)
           {print(c('warning at iter:', s)); print(war);traceback()},
        error =function(e)
            {print(c('error at iter:', s)); print(e);traceback()} ,
        finally = print("Done!"))

    if(burn.in==0) burn.in = 1 else burn.in = 1:(burn.in+1)
    y0 =  y0[,-burn.in]
    w0 =  w0[,-burn.in]
    if(uncertainty)  g0 = g0[-burn.in] else g0 = NULL
    list(w = w0, y = y0, g = g0, burn.in = max(burn.in)-1,
         sd = list(w=w_sd, y = y_sd))
}

### ==================================================
### Supporting functions
### ==================================================
rExp<-function(l,a=1){
    ## Sampling from a truncated exponential distribution
 	unif = runif(length(l))
	-log(1- unif*(1-exp(-l*a)))/(l + tol.err)
}

rExp2<-function(l, g, Z, Z0){
    ## Sampling from a zero-inflated Gumbel or a one-inflated Exponential
    ## Method one: P(z=0|g) = 1\delta_{s=0} + g\delta_{s>0}; dirac delta
    p = 1- exp(-l)
    unif = matrix(runif(length(l)), dim(p))
    U = 1 + 0*p
    U[!Z0] = -log(1- unif[!Z0]*p[!Z0])/(l[!Z0] + tol.err)
    aa = (unif < g*p/(g*p + 1-p)) & Z0
    U[aa] =  -log(1 - unif[aa]*(g*p[aa] + 1-p[aa])/g)/l[aa]
    U
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
    eta.prop = eta.old + rnorm(length(eta.old), 0, sd = eta_sd)
    dist = 1/cophenetic(eb.phylo(tree, tree.ht, eta.prop))
    diag(dist)<-0
    pdist.new = c(dist[i,]%*%Z)
    no0 = which(pd0)
    if(length(no0)){
        if(length(no0)==sum(Z[i,])) likeli = -Inf else 
        likeli = sum(((log(pdist.new)- log(pdist.old))*Z[i,])[-no0])-
            sum((ywU*(pdist.new - pdist.old))[-no0]) 
    }else{
        likeli = sum((log(pdist.new)- log(pdist.old))*Z[i,] )-
            sum(ywU*(pdist.new - pdist.old))
    }
    ratio = min(1, exp(likeli));ratio
    u = (runif(1)<=ratio)
    if(u) { eta.old  = eta.prop; pdist.old = pdist.new}
    
    list (eta=eta.old, dist=pdist.old)
}

rg<-function(Z,l){
    ## Sampling the uncertainty parameter
    ## Method one: P(z=0|g) = 1\delta_{s=0} + g\delta_{s>0}; dirac delta
    ZZ = 1*(l<1)                  # 1-inflated exponential
    ## Method 1
    M = sum(ZZ*Z)              # l<1(OR S>0) and Z=1,  N++
    N = sum((1-Z)*ZZ)          # l<1(OR S>0) and Z=0,  N-+
    rbeta(1 , N + 1, M + 1)         
}


lof<-function(Z, indeces = FALSE){
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
    if(indeces) bank else  Z[,bank]
}

### ==================================================
### Tree functions
### ==================================================
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

eb.phylo<-function(phy, heights, a){
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

heights.phylo<-function(x){
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

### ==================================================
