

ICM_est_C<-function(Z, tree, slices = 10, distOnly = FALSE, uncertainty = FALSE, sparse=FALSE, ...){
    ## loading extra args
    el <-list(...)
    ## parameters set-up
    nw = ncol(Z);ny = nrow(Z)
    n = nw*ny                          

    y = if(!is.null(el$y)) el$y else 1
    w = if(!is.null(el$w)) el$w else 1
    a_w = if(!is.null(el$a_w)) el$a_w else 1
    a_y = if(!is.null(el$a_w)) el$a_y else  1
    b_w = b_y = 1;
    y_sd = if(is.null(el$y_sd)) 0.2 else el$y_sd
    w_sd = if(is.null(el$w_sd)) 0.2 else el$w_sd
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
    w0in = matrix(w, nrow= nw, ncol =subItra+1)
    y0in = matrix(y, nrow= ny, ncol =subItra+1)
    petain = rep(0, subItra+1)
    g0in = rep(0, subItra+1)

    Z0 <- which(Z==0)
	mc <- colSums(Z)
	mr <- rowSums(Z)

    lowerIndex = lower.tri.Index(ny)
    upperIndex = upper.tri.Index(ny)

    ## Arranging the tree
    tree.ht = arrange.tree(tree)
    tree$tip.label = 1:length(tree$tip.label)
    dist = cophFast(eb.phylo(tree, tree.ht, peta[1]), lowerIndex, upperIndex, ny)
    sparseZ = Z
    if(sparse){
        ind <- which(Z==1, arr.ind=TRUE)
        sparseZ = sparseMatrix(ind[,1], j=ind[,2], x= rep(1, nrow(ind)),  dims=dim(Z))
    }
    ind <- seq.int(1, prod(dim(sparseZ)),  nrow(sparseZ))
    pdist = dist%*%sparseZ
    if(sparse)  pdist = as(pdist, 'matrix')
    pdist0 = apply(pdist, 1, function(r) which(r==0))
    pdist00 = which(pdist==0)
    
    print(sprintf("Run for %i slices with %i burn-ins",slices, burn.in))

    i=1;s=1
    tryCatch(
        for(s in 1:slices){
            if(s%%100==0)
                print(sprintf('slice: %d', s))
            ## Arranging the tree
            dist =cophFast(eb.phylo(tree, tree.ht, peta[s]),lowerIndex, upperIndex, ny)
            pdist = dist%*%sparseZ
            if(sparse)
                pdist = as(pdist, 'matrix')
            pdist[pdist00]<-1
            
            ## Updating latent scores
            yw = outer(y0[,s],w0[, s])
            if(!uncertainty){
                U0 <- rExp(pdist*yw)
                U0[Z0] <- 1
            }else
                U0 <-rExp2(pdist*yw, g0[s], Z, Z0)
        
            ## looping over nrow(Z), for ICM
            Upd = U0*pdist
            for (i in 1:subItra){
                if(!distOnly){                
                    ## Updating the parasite parameters
                    w0in[, i+1]<-raffinity.MH(w0in[,i],Z[i,],
                                              y0in[i,i]*Upd[i,],
                                              sig=w_sd, c(a_w, 1))
                    ## Updating host parameters
                    y0in[, i+1]<-raffinity.MH(y0in[,i],mr,
                                              tcrossprod(w0in[,i+1], Upd),
                                              sig=y_sd, c(a_y, 1))
                }
                ## Updating similarity matix parameter
                new.eta = rEta.copheneticFastSpa(petain[i],tree,tree.ht,
                    pdist[i,],pdist0[[i]],i,sparseZ,Z,
                    y0in[i,i+1]*(w0[,s]*U0[i,]),
                    eta_sd, lowerIndex, upperIndex, ny, ind, sparse)
                petain[i+1] = new.eta$eta
                if(new.eta$change)
                    pdist[i,] = new.eta$dist # very slow when using sparse assignment
                
                if(uncertainty){
                    g0in[i+1] = rg(Z[i,], l=U0[i,])
                }
            }
            
            ## MH Adaptiveness
            ## Parasite parameters (w)
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


rEta.copheneticFastSpa<-function(eta.old,tree,tree.ht,pdist.old, no0,i, sZ, Z, ywU,
                              eta_sd =0.01, a, b,nr, ind, sparse){
    ## a faster version of rEta using faster cophenetic and sparse matrices
    change= FALSE
    eta.prop = eta.old + eta_sd*rnorm(1)
    dist = cophFast(eb.phylo(tree, tree.ht, eta.prop), a, b,nr)
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
    if(runif(1)<= min(1, exp(likeli)))
        { eta.old  = eta.prop; pdist.old = c(pdist.new);change=TRUE}
    
    list (eta=eta.old, dist=pdist.old, change=change)
}

