#!/bin/R
## Tolerance error 
tol.err=1e-5

### ==================================================
### Main MCMC function
### ==================================================
network_est<-function(Z, slices = 10, tree = NULL, model.type = c('full', 'distance', 'affinity'), uncertainty = FALSE,sparse=TRUE, ... ){
    require(geiger)
    require(phangorn)
    require(Matrix)
    require(methods)
    ## Running options:
    ## 1 - Full to combined 2 and 3;
    ## 2 - Affinity-only model (affinity);
    ## 3 - distance-only model (distance).
    model.type= tolower(model.type[1])
    ## General warnings are checks
    if(missing(Z)) stop('Interaction matrix is missing!')
    if(slices==0)  stop('no. of slices cannot be 0!')

    cleaned = network_clean(Z, tree, model.type)
    Z = cleaned$Z
    tree = cleaned$tree

    ## For affinity-only model 
    if(grepl('aff', model.type)){
        ## sparse option is not used for affinity
        print('Running affinity model...')
        param = fullJoint_est(Z, iter = slices, uncertainty = uncertainty, ...)
    }
    ##  Full and distance model
    if(grepl('(dist|full)', model.type)){
        ## Running the MCMC
        print(paste0('Running ',
                     ifelse(grepl('dist', model.type),
                            'distance model...', 'full model...')))
        param  = ICM_est(unname(Z),tree,slices, distOnly = grepl('dist', model.type),
            uncertainty = uncertainty, sparse=sparse, ...)
    }
    list(param =param , Z = Z, tree=tree)
}

network_clean<-function(Z, tree = NULL, model.type = c('full', 'distance', 'affinity'), uncertainty = FALSE){
    ## A function to clean Z and tree as:
    ## - converts Z to binary, if uncertain is FALSE
    ## - removed empty columns from Z
    ## - confirms tree is a phylo object when full or distance models
    ## - removed tree tips that do not exist in the rows of Z
    ## - removes the rows of Z with no corresponding tree tips
    ## - orders the rows of Z abiding to cophenetic conversion.
    
    require(geiger)
    require(phangorn)
    require(Matrix)
    ## Running options:
    ## 1 - Full to combined 2 and 3;
    ## 2 - Affinity-only model (affinity);
    ## 3 - distance-only model (distance).

    model.type= tolower(model.type[1])
    ## General warnings are checks
    if(missing(Z)) stop('Interaction matrix is missing!')
    if(!all(range(Z)==c(0,1)) & !uncertainty){
        warning('Z is converted to binary!', immediate. = TRUE, call.= FALSE)
        Z = 1*(Z>0)
    }
    if(any(colSums(Z)==0)){
        stop('Z has empty columns, please remove!', immediate. = TRUE, call.= FALSE)
        Z = Z[,which(colSums(Z)>0)]
    }
    if(grepl('aff', model.type)){
        if(!is.null(tree))
            warning('affinity model is chosen; ignoring tree!',
                    immediate. = TRUE, call. = FALSE)
    }
    if(grepl('(dist|full)', model.type)){
        if(is.null(tree))
            stop('distance model is chosen, but tree is null!')
        ## testing tree is a phylo object
        if(!is.phylo(tree))
            stop('tree must be a phylogeny tree, see gieger!')

        ## testing that all tips exist in Z
        if(!all(tree$tip.label %in% rownames(Z))){
            warning('not all species in tree exist in Z; missing are removed from tree!',
                    immediate.= TRUE, call. = FALSE)
            tree = drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(Z)])
        }
        ## Testing all names in com exist in dist
        if(!all(rownames(Z) %in% tree$tip.label)){
            warning('not all row-species in Z exist tree; missing are removed from Z!',
                    immediate.= TRUE, call. = FALSE)
            Z = Z[rownames(Z) %in% tree$tip.label,]
        }
        ## Testing that names in Z exist only once
        if(!all(sapply(sapply(rownames(Z),
                              function(r) which(r==tree$tip.label)), length)==1))
            stop('some row-species in Z exist more than once tree!')
        
        ## Making sure the order of hosts in Z and tree are the same
        aux  = cophenetic(tree)
        if(!all(rownames(aux)==rownames(Z))){
            row.order <- sapply(rownames(aux),function(r) which(r==rownames(Z)))
            print('Ordering the rows of Z to match tree...')
            Z = Z[row.order,]
        }
        if(max(range(tree$edge.length))>1){
            print('normalizing tree edges by the maximum pairwise distance!')
            tree$edge.length = tree$edge.length/max(aux)
        }
    }
    list(Z = Z, tree = tree)
}

ICM_est<-function(Z, tree, slices = 10, distOnly = FALSE, uncertainty = FALSE, sparse=TRUE, ...){
    ## Main function for ICM, only applied is when distance input is used
    ## loading extra args
    el <-list(...)
    ## parameters set-up
    nw = ncol(Z);ny = nrow(Z)
    n = nw*ny                          

    y = if(!is.null(el$y)) el$y else 1
    w = if(!is.null(el$w)) el$w else 1
    a_w = if(!is.null(el$a_w)) el$a_w else 1
    a_y = if(!is.null(el$a_y)) el$a_y else  1
    b_w = if(!is.null(el$b_w)) el$b_w else 1
    b_y = if(!is.null(el$b_y)) el$b_y else 1
    y_sd = if(is.null(el$y_sd)) rep(0.2, ny) else el$y_sd
    w_sd = if(is.null(el$w_sd)) rep(0.2, nw) else el$w_sd
    eta = if(is.null(el$eta)) 0 else el$eta
    eta_sd = if(!is.null(el$eta_sd)) el$eta_sd else 0.005

    ## Burn in set-up
    beta = if(!is.null(el$beta)) el$beta else 1
    burn.in = if(is.null(el$burn.in)) floor(0.5*slices) else floor(el$burn.in*slices)

    ## Variable holders
    ## outer loop
    y0<-matrix(y, nrow= ny, ncol =slices+1)
    w0<-matrix(w, nrow= nw, ncol =slices+1)
    g0<-rep(0, slices+1)
    peta = rep(eta, slices)

    ## inner loop
    subItra = ny 
    g0in = rep(0, subItra+1)

    ## New code to same memory
    w0.new = w0.count = w0.sum=0
    y0.new = y0.count = y0.sum=0
    w0.last = if(length(w)==1) rep(w, nw) else w
    y0.last = if(length(y)==1) rep(y, ny) else y
    peta.new = peta.count = peta.sum =0
    peta.last = eta
    
    ## special variables
    Z00 = Z==0
    Z0 <- which(Z==0)
	mc <- colSums(Z)
	mr <- rowSums(Z)

    ## indexing the upper and lower triangular of the matrix
    lowerIndex = lower.tri.Index(ny)
    upperIndex = upper.tri.Index(ny)

    ## Arranging the tree
    tree.ht = arrange.tree(tree)
    tree$tip.label = 1:length(tree$tip.label) # removing tip labels

    ##    dist = cophFast(eb.phylo(tree, tree.ht, peta[1]), lowerIndex, upperIndex, ny)
    dist = 1/cophenetic(eb.phylo(tree, tree.ht, peta[1]))
    diag(dist)<-0
    
    sparseZ = Z
    if(sparse){
        ind <- which(Z==1, arr.ind=TRUE)
        sparseZ = sparseMatrix(ind[,1], j=ind[,2], x= rep(1, nrow(ind)),  dims=dim(Z))
    }
    ind <- seq.int(1, prod(dim(sparseZ)),  nrow(sparseZ))
    pdist = dist%*%sparseZ
    if(sparse)  pdist = as(pdist, 'matrix')
    pdist0 = apply(pdist, 1, function(r) which(r==0))
    if(length(pdist0) && length(pdist0)!=ny) stop('pdist0 error')
    pdist00 = which(pdist==0)

    ## starting the loop
    print(sprintf("Run for %i slices with %i burn-ins",slices, burn.in))
    print(paste('Matrix dimension:', nrow(Z)[1],'x', ncol(Z) ))
    i=1;s=1
    tryCatch(
        for(s in 1:slices){
            if(s%%200==0){
                print(sprintf('slice: %d, at %s', s, Sys.time())) # iteration update
                if(!is.null(el$backup))
                    save(y0,w0,peta,g0, file='snapshot.RData')
            }
            
            ## Updating the tee
            ## dist =cophFast(eb.phylo(tree, tree.ht, peta[s]),lowerIndex, upperIndex, ny)
            dist = 1/cophenetic(eb.phylo(tree, tree.ht, peta[s]))
            diag(dist)<-0
            pdist = dist%*%sparseZ
            ## conversion takes less time then extraction/insertion from dgMatrix class
            if(sparse) pdist = as(pdist, 'matrix') 
            pdist[pdist00]<-1
            
            ## Updating latent scores
            yw = outer(y0[,s],w0[, s])
            if(!uncertainty){
                U0 <- rExp(pdist*yw)
                U0[Z0] <- 1
            }else
                U0 <-rExp2(pdist*yw, g0[s], Z,Z0, Z00)

            ## looping over nrow(Z), for ICM
            Upd = U0*pdist
            for (i in 1:subItra){
                if(!distOnly){                
                    ## Updating the parasite parameters
                    w0.new<-raffinity.MH(w0.last,Z[i,],
                                         y0.last[i]*(Upd[i,]),
                                         sig=w_sd, c(a_w, b_w))
                    w0.count = w0.count + 1*(abs(w0.new-w0.last) > tol.err)
                    w0.sum = w0.sum + w0.new*mr[i]
                    w0.last = w0.new
                    ## Updating host parameters
                    y0.new<-raffinity.MH(y0.last,mr,
                                         tcrossprod(w0.new,Upd),
                                         sig=y_sd, c(a_y, b_y))
                    y0.count = y0.count + 1*(abs(y0.new-y0.last) > tol.err)
                    y0.sum = y0.sum + y0.new
                    y0.last = y0.new
                    
                }
                ## Updating similarity matix parameter
                new.eta = rEta.copheneticFast(peta.last,tree,tree.ht,
                    pdist[i,],
                    if(length(pdist0)) pdist0[[i]] else NULL,i,sparseZ,Z,
                    y0.new[i]*(w0[,s]*U0[i,]),
                    eta_sd, lowerIndex, upperIndex, ny, ind, sparse)
                peta.new = new.eta$eta
                peta.count = peta.count + 1*(abs(peta.new - peta.last) > tol.err)
                peta.sum = peta.sum + peta.new*mr[i]
                peta.last = peta.new
                if(new.eta$change)
                    pdist[i,] = new.eta$dist # very slow when using sparse assignment
                
                if(uncertainty){
                    g0in[i+1] = rg(Z[i,], l=U0[i,])
                }
            }
            
            ## MH Adaptiveness
            # Parasite parameters (w)
            w_sd = w_sd*exp(beta*(w0.count/subItra-0.44)/log(1 +s))
            y_sd = y_sd*exp(beta*(y0.count/subItra-0.44)/log(1 +s))
            ## Tree scaling parameter (eta)
            eta_sd = eta_sd*exp(beta*(peta.count/subItra-0.44)/log(1 +s))
            
            if(!distOnly){
                w0[,s+1] <- w0.sum/sum(mr)
                w0.new = w0.count = w0.sum=0
                w0.last = w0[,s+1]

                y0[,s+1]<- y0.sum/subItra
                y0.new = y0.count = y0.sum=0
                y0.last = y0[,s+1]
            }
            peta[s+1]<- peta.sum/sum(mr)
            peta.new = peta.count = peta.sum =0
            peta.last = peta[s+1]

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

    y = if(!is.null(el$y)) el$y else 1
    w = if(!is.null(el$w)) el$w else 1
    a_w = if(!is.null(el$a_w)) el$a_w else 1
    a_y = if(!is.null(el$a_y)) el$a_y else  1
    b_w = if(!is.null(el$b_w)) el$b_w else 1
    b_y = if(!is.null(el$b_y)) el$b_y else 1
    y_sd = if(is.null(el$y_sd)) 0.2 else el$y_sd
    w_sd = if(is.null(el$w_sd)) 0.2 else el$w_sd
    
    ## Burn in set-up
    batch.size = if(!is.null(el$batch.size)) el$batch.size else 50
    beta = if(!is.null(el$beta)) el$beta else 1
    burn.in = if(is.null(el$burn.in)) floor(0.5*iter) else floor(el$burn.in*iter)

    ## Variable holders
    y0<-matrix(y, nrow= ny, ncol =iter+1)
    w0<-matrix(w, nrow= nw, ncol =iter+1)
    g0<-rep(0, iter+1)
    peta = rep(0, iter+1)
    
    Z0 <-Z==0
	mc <- colSums(Z)
	mr <- rowSums(Z)

    print(sprintf("Run for %i iterations with %i burn-ins",iter, burn.in))
    s=1
    tryCatch(
        for(s in 1:iter){
            if(s%%200==0){
                print(sprintf('iteration %d, at %s', s, Sys.time()))
                if(!is.null(el$backup))
                    save(y0,w0,peta,g0, file='snapshot.RData')
            }
            ## Updatting latent scores
            U0 <- rExp(outer(y0[,s],w0[, s]))
            U0[Z0] <- 1
            
            ## Updating the parasite parameters
            w0[, s+1]<-raffinity.MH(w0[,s],mc,
                                    crossprod(y0[,s],U0),
                                    sig=w_sd, c(a_w, b_w))
            ## Updating host parameters
            y0[, s+1]<-raffinity.MH(y0[,s],mr,
                                    tcrossprod(w0[,s+1],U0),
                                    sig=y_sd, c(a_y, b_y))
            ## Uncertain parameter sampling
            
            
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

rExp2<-function(l, g, Z, Z0, Z00){
    ## Sampling from a zero-inflated Gumbel or a one-inflated Exponential
    ## Method one: P(z=0|g) = 1\delta_{s=0} + g\delta_{s>0}; dirac delta
    p = 1- exp(-l)
    unif = matrix(runif(length(l)), dim(p))
    U = 1 + 0*p
    U[-Z0] = -log(1- unif[-Z0]*p[-Z0])/(l[-Z0] + tol.err)
    aa = Z00 & (unif < g*p/(tol.err + g*p + 1-p))
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
    pdist.new = c(crossprod(dist[i,],Z))
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

### ==================================================
### Tree functions
### ==================================================
arrange.tree<-function(phy){
    ## taken from .b.phylo from
    ## https://github.com/mwpennell/geiger-v2/blob/master/R/utilities-phylo.R#L1353-L1382
    ht=heights.phylo(phy)
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
### Faster Cophenetic using Phangorn package and ape>5.0
### ==================================================

lower.tri.Index <-function(n){
    ## returning the vectorized indices of the lower triangular
    ## elements of the matrix, excluding diagonal 
    unlist(sapply(1:(n-1), function(i) n*(i-1) + (i+1):n))
}

upper.tri.Index<-function(n){
    ## returning the vectorized indices of the upper triangular
    ## elements of the matrix, excluding diagonal 
    b = matrix(0, n ,n )
    lowerIndex = lower.tri.Index(n)
    b[lowerIndex]<-lowerIndex
    b =t(b)
    a = which(upper.tri(b))
    a[order(b[upper.tri(b)])]
}


cophFast<-function(tree, lowerIndex, upperIndex, n){
    ## a fast version of phylo tree using phangorn:::coph
    ## phangorn:::coph returns a dist object.
    b = matrix(0, n,n)
    b[lowerIndex]<- b[upperIndex]<-(1/phangorn:::coph(tree))
    b
}

rEta.copheneticFast<-function(eta.old,tree,tree.ht,pdist.old, no0,i, sZ, Z, ywU,
                              eta_sd =0.01, a, b,nr, ind, sparse){
    ## a faster version of rEta using faster cophenetic and sparse matrices
    change = FALSE
    eta.prop = eta.old + eta_sd*rnorm(1)
    dist = cophFast(eb.phylo(tree, tree.ht, eta.prop), a, b,nr)
    if(any(is.na(dist))){
        return( list (eta=eta.old, dist=pdist.old, change=change))
    }else{
    pdist.new = dist[i,]%*%sZ
    if(sparse)
        pdist.new= pdist.new@x
    if(length(no0)){
        if(length(no0)==sum(Z[i,])) likeli = -Inf else 
        likeli = sum(((log(pdist.new)- log(pdist.old))*Z[i,])[-no0])-
            sum((ywU*(pdist.new - pdist.old))[-no0]) 
    }else{
        likeli = sum((log(pdist.new)- log(pdist.old))*Z[i,] )-
            sum(ywU*(pdist.new - pdist.old))
    }
    if(!is.na(likeli) && !is.nan(likeli) && runif(1)<= min(1, exp(likeli)))
        { eta.old  = eta.prop; pdist.old = c(pdist.new);change=TRUE}
    
    list (eta=eta.old, dist=pdist.old, change=change)
}
}

### ==================================================
### Faster extraction for Matrix package, class dgCMatrix
### ==================================================

row_extract<-function(m, i, ind){
    ## extracting row i from matrix m based on ind
    ## m: matrix of package Matrix
    ## i: is the row number
    ## ind: is the index of the first cell of each column
    ## ind can be constructed as
    ## seq.int(1, prod(dim(m)),  nrow(m))
    
    if(class(m)=='dgeMatrix')
        m@x[ind + (i-1)] else m[ind + (i-1)]
}


### ==================================================
