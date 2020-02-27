ICM_est_over_acc <-
function(Z, tree, slices = 10, distOnly = FALSE, uncertainty = FALSE, sparse=TRUE, ...){
    ## Main function for ICM, only applied is when distance input is used
    ## loading extra args
    el <-list(...)
    ## parameters set-up
    nw = ncol(Z);ny = nrow(Z)
    n = nw*ny                          
    tol.err = 1e-4
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
    batch.size = if(!is.null(el$batch.size)) el$batch.size else 50
    beta = if(!is.null(el$beta)) el$beta else 1
    burn.in = if(is.null(el$burn.in)) floor(0.5*slices) else floor(el$burn.in*slices)

    ## Variable holders
    ## outer loop
    y0<-matrix(y, nrow= ny, ncol =slices+1)
    w0<-matrix(w, nrow= nw, ncol =slices+1)
    g0<-rep(0, slices+1)
    peta = rep(eta, slices)

    ## inner loop
    ## g0in = rep(0, subItra+1)

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
    t.max = get.max.depth(tree)
    ## tree$tip.label = 1:length(tree$tip.label) # removing tip labels
    dist.original = unname(cophenetic(rescale(tree, 'EB', 0)))/2
    dist = 1/EB.distance(dist.original, t.max, peta[1])
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
    tryCatch(
        for(s in 1:slices){
            if(s%%200==0)
                print(sprintf('slice: %d, at %s', s, Sys.time())) # iteration update

            ## Updating the tee
            dist = 1/EB.distance(dist.original, t.max, peta[s])
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
            Upd = U0 *pdist
            if(!distOnly){
                ## Updating the parasite parameters
                w0[, s+1]<-raffinity.MH(w0[,s],mc,
                                        crossprod(y0[,s],U0),
                                        sig=w_sd, c(a_w, b_w))
                ## Updating host parameters
                y0[, s+1]<-raffinity.MH(y0[,s],mr,
                                        tcrossprod(w0[,s+1],U0),
                                        sig=y_sd, c(a_y, b_y))
                ## Uncertain parameter sampling
            }
            ## Updating similarity matix parameter
            new.eta = rEta.over.acc(peta[s],pdist,pdist0,pdist00, sparseZ,
                Z,outer(y0[,s+1], w0[,s+1])*U0,eta_sd,sparse,
                                    dist.original, t.max)
            peta[s+1] = new.eta$eta
            if(new.eta$change)
                pdist = new.eta$dist # very slow when using sparse assignment
            
            if(uncertainty)
                g0[s+1] = rg(Z, l=U0)
            
            ## MH Adaptiveness
            ## Parasite parameters (w)
            if(s%%batch.size==0){
                ss = s+1 - 1:batch.size
                ac =1- rowMeans(abs(w0[,ss] - w0[,ss+1])<tol.err)
                w_sd = w_sd*exp(beta*(ac-0.44)/log(1 +s/batch.size))
                
                ac =1- rowMeans(abs(y0[,ss] - y0[,ss+1])<tol.err)
                y_sd = y_sd*exp(beta*(ac-0.44)/log(1 +s/batch.size))

                ac =1- mean(abs(peta[ss] - peta[ss+1])<tol.err)
                eta_sd = eta_sd*exp(beta*(ac-0.44)/log(1 +s/batch.size))
            }
            
            ## Tree scaling parameter (eta)
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
