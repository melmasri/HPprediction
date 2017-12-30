library(phangorn)
library(Matrix)


ICM_est_C<-function(Z, tree, slices = 10, distOnly = FALSE, uncertainty = FALSE,sparse=FALSE, ...){
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
    y_sd = if(is.null(el$y_sd)) rep(0.2, ny) else el$y_sd
    w_sd = if(is.null(el$w_sd)) rep(0.2, nw) else el$w_sd
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
    dist = 1/cophenetic(eb.phylo(tree, tree.ht, peta[1]))
    diag(dist)<-0
    sparseZ = Z
    if(sparse){
        ind <- which(Z==1, arr.ind=TRUE)
        sparseZ = sparseMatrix(ind[,1], j=ind[,2], x= rep(1, nrow(ind)),  dims=dim(Z))
    }
    ind <- seq.int(1, prod(dim(sparseZ)),  nrow(sparseZ))
    pdist = dist%*%Z
    ## pdist0 = pdist==0
    pdist = apply(pdist, 1, function(r) which(r==0))
    pdist00 = which(pdist0)
    
    lowerIndex = lower.tri.Index(ny)
    upperIndex = upper.tri.Index(ny)

    U0=matrix(1,ny, nw)
    i=1;s=1
    tryCatch(
        for(s in 1:slices){
            if(s%%100==0)
                print(sprintf('slice: %d', s))
            ## Arranging the tree
            dist = 1/cophFast(eb.phylo(tree, tree.ht, peta[s]),
                lowerIndex, upperIndex, ny)
            diag(dist)<-0
            pdist = dist%*%sparseZ
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
                                              y0in[i,i]*row_extract(Upd, i, ind),
                                              sig=w_sd, c(a_w, 1))
                    ## Updating host parameters
                    y0in[, i+1]<-raffinity.MH(y0in[,i],mr,
                                              tcrossprod(w0in[,i+1], Upd),
                                              sig=y_sd, c(a_y, 1))
                }
                ## Updating similarity matix parameter
                new.eta = rEta.copheneticFast(petain[i],tree,tree.ht,
                    row_extract(pdist,i,ind),pdist0[[i]],i,sparseZ,Z,
                    y0in[i,i+1]*(w0[,s]*row_extract(U0,i, ind)),
                    eta_sd, lowerIndex, upperIndex, ny, ind)
                petain[i+1] = new.eta$eta
                pdist[i,] <-new.eta$dist
                
                if(uncertainty){
                    g0in[i+1] = rg(Z[i,], l=row_extract(U0,i, ind))
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


row_extract<-function(m, i, ind){
    ## extracting row i from matrix m based on ind
    if(class(m)=='dgeMatrix')
        m@x[ind + (i-1)] else m[ind + (i-1)]
}

rEta.copheneticFast<-function(eta.old,tree,tree.ht,pdist.old, pd0,i, sZ, Z, ywU,
                              eta_sd =0.01, a, b,nr, ind){
    ## a faster version of rEta using faster cophenetic and sparse matrices
    eta.prop = eta.old + eta_sd*rnorm(1)
    ##dist = 1/cophenetic(eb.phylo(tree, tree.ht, eta.prop))
    dist = 1/cophFast(eb.phylo(tree, tree.ht, eta.prop), a, b,nr)
    diag(dist)<-0
    pdist.new = dist[i,]%*%sZ
    ## no0 = which(pd0)
    if(length(no0)){
        if(length(no0)==sum(Z[i,])) likeli = -Inf else 
        likeli = sum(((log(pdist.new)- log(pdist.old))*Z[i,])[-no0])-
            sum((ywU*(pdist.new - pdist.old))[-no0]) 
    }else{
        likeli = sum((log(pdist.new)- log(pdist.old))*Z[i,] )-
            sum(ywU*(pdist.new - pdist.old))
    }
    ## ratio = min(1, exp(likeli))
    u = runif(1)<= min(1, exp(likeli))
    if(u) { eta.old  = eta.prop; pdist.old = pdist.new}
    
    list (eta=eta.old, dist=pdist.old)
}

