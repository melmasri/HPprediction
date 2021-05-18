ICM_est <-
function(Z, tree, slices = 10, distOnly = FALSE, uncertainty = FALSE, sparse=TRUE, ...){
    ## Main function for ICM, only applied is when distance input is used
    ## loading extra args
    el <-list(...)
    ## parameters set-up
    nw = ncol(Z);ny = nrow(Z)
    n = nw*ny                          
    tol.err = 1e-8

    meta = distances_metadata(tree)
    tree = meta$dist
    kernel_name = meta$kernel
    kernel_func =meta$kernel_func
    t.max = meta$t.max


    y = if(!is.null(el$y)) el$y else 1
    w = if(!is.null(el$w)) el$w else 1
    a_w = if(!is.null(el$a_w)) el$a_w else 1
    a_y = if(!is.null(el$a_y)) el$a_y else  1
    b_w = if(!is.null(el$b_w)) el$b_w else 1
    b_y = if(!is.null(el$b_y)) el$b_y else 1
    y_sd = if(is.null(el$y_sd)) rep(0.2, ny) else el$y_sd
    w_sd = if(is.null(el$w_sd)) rep(0.2, nw) else el$w_sd
    eta = if(is.null(el$eta)) 0.1 else el$eta
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
    N = outer(rowSums(Z), colSums(Z))
    
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
    dummyMat <- matrix(0, ny, ny)
    ## Arranging the tree
    ## tree.ht = arrange.tree(tree)
    if(is.phylo(tree))
        dist.original = unname(cophenetic(rescale(tree, 'EB', 0)))/2
    else
        dist.original = tree
    ##dist = 1/kernel_func(dist = dist.original, tmax = t.max, eta = peta[1])
    dist = kernel_func(dist = dist.original, tmax = t.max, eta = peta[1])
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
            dist = kernel_func(dist = dist.original, tmax = t.max, eta=peta[s])
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
                    w0.sum = w0.sum + w0.new
                    w0.last = w0.new
                    ## Updating host parameters
                    y0.new<-raffinity.MH(y0.last,mr,
                                         tcrossprod(w0.new,Upd),
                                         sig=y_sd, c(a_y, b_y))
                    y0.count = y0.count + 1*(abs(y0.new-y0.last) > tol.err)
                    y0.sum = y0.sum + y0.new
                    y0.last = y0.new
                    
                }## Updating similarity matix parameter
                new.eta = rEta(peta.last,
                               pdist[i,],
                               if(length(pdist0)) pdist0[[i]] else NULL,i,sparseZ,Z,
                               if(distOnly) U0[i,]else  y0.new[i]*(w0.new*U0[i,]),
                               eta_sd,
                               sparse,
                               dist.original,
                               t.max,
                               kernel_func,
                               kernel_name)
                peta.new = new.eta$eta
                peta.count = peta.count + 1*(abs(peta.new - peta.last) > tol.err)
                peta.sum = peta.sum + peta.new
                peta.last = peta.new
                ## if(new.eta$change){
                ##     pdist[i,] = new.eta$dist # very slow when using sparse assignment
                ## }
                
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
                w0[,s+1] <- w0.sum/subItra
                w0.new = w0.count = w0.sum=0
                w0.last = w0[,s+1]

                y0[,s+1]<- y0.sum/subItra
                y0.new = y0.count = y0.sum=0
                y0.last = y0[,s+1]
            }
            peta[s+1]<- peta.sum/subItra
            peta.new = peta.count = peta.sum =0
            peta.last = peta[s+1]

            if(uncertainty){
                g0[s+1] = sum(g0in[-1])/subItra
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
         sd = list(w= if(distOnly) NULL else w_sd, y = if(distOnly) NULL else y_sd, eta= eta_sd))
}
