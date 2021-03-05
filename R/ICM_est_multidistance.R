

ICM_est_multidistance <-function(Z,
                                 distances,
                                 slices = 10,
                                 distOnly = FALSE,
                                 uncertainty = FALSE,
                                 sparse=TRUE,
                                 ...){
    ## Main function for ICM, only applied is when distance input is used
    ## loading extra args
    el <-list(...)
    
    ## number of distances
    num.dist = length(distances)
    
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
    eta = if(is.null(el$eta)) rep(0, num.dist) else el$eta
    eta_sd = if(!is.null(el$eta_sd)) el$eta_sd else 0.005
    dist.weights  = if(!is.null(el$dist.weights)) el$dist.weights else  rep(1, num.dist)/num.dist
    weights_sd = if(is.null(el$weights_sd)) 0.2 else el$weights_sd[1]
    weights_hyper = if(is.null(el$weights_hyper)) 10*rep(1, num.dist)/num.dist else sum(el$weights_hyper)*rep(1, num.dist)/num.dist

    if(length(eta)!=num.dist) eta = rep(eta[1], num.dist)
    if(length(eta_sd) != num.dist) eta_sd = rep(eta_sd[1], num.dist)
    if(length(dist.weights)!=num.dist || sum(dist.weights)!=1){
        warning('dist.weights either do not sum to 1, or are not equal in length to the number of passed distances, converting to uniform weights',
                immediate.= TRUE, call. = FALSE)
        dist.weights = rep(1, num.dist)/num.dist
    }

    ## Burn in set-up
    beta = if(!is.null(el$beta)) el$beta else 1
    burn.in = if(is.null(el$burn.in)) floor(0.5*slices) else floor(el$burn.in*slices)
    
    ## Variable holders
    ## outer loop
    y0<-matrix(y, nrow= ny, ncol =slices+1)
    w0<-matrix(w, nrow= nw, ncol =slices+1)
    g0<-rep(0, slices+1)
    peta = matrix(eta, nrow = num.dist, ncol = slices + 1)
    
    ## inner loop
    subItra = ny 
    g0in = rep(0, subItra+1)

    ## New code to same memory
    w0.new = w0.count = w0.sum=0
    y0.new = y0.count = y0.sum=0
    w0.last = if(length(w)==1) rep(w, nw) else w
    y0.last = if(length(y)==1) rep(y, ny) else y
    peta.new = peta.count = peta.sum =numeric(num.dist)
    peta.last = eta

    dist.weights.count =  dist.weights.sum = 0
    dist.weights.new = dist.weights.last = dist.weights
    dist.weights = matrix(0, nrow = num.dist, ncol = slices +1)
    dist.weights[,1] <-dist.weights.new

    ## indexing the upper and lower triangular of the matrix
    lowerIndex = lower.tri.Index(ny)
    upperIndex = upper.tri.Index(ny)
    dummyMat <- matrix(0, ny, ny)
    ## Arranging the tree
    ## tree.ht = arrange.tree(tree)

    t.max = lapply(distances, get.max.depth)
    ## tree$tip.label = 1:length(tree$tip.label) # removing tip labels
    
    dist.original = lapply(distances,
                           function(r)
                               if(is.phylo(r))
                                   unname(cophenetic(rescale(r, 'EB', 0))/2)
                               else
                                   unname(r)
                           )
    
    ## special variables
    Z00 = Z==0
    Z0 <- which(Z==0)
	mc <- colSums(Z)
	mr <- rowSums(Z)

    transform_dist <-function(d, tm, eta){
        aa = lapply(1:num.dist, function(i){
            a = 1/(EB.distance(d[[i]], tm[[i]], eta[i]))
            diag(a) <-0
            a
        })
        aa
    }
    
    transform_dist.inv <-function(d, tm, eta){
        aa = lapply(1:num.dist, function(i){
            EB.distance(d[[i]], tm[[i]], eta[i])
        })
        aa
    }
    
    
    ## dist.list.org = transform_dist(dist.original, t.max, peta[,1])
    ## dist = matrix(matrix(unlist(dist.list.org),
    ##                      ny*ny, num.dist) %*% dist.weights.last, ny, ny)

    dist.list.org.inv = transform_dist.inv(dist.original, t.max, peta[,1])

    dist.inv = 1/matrix(matrix(unlist(dist.list.org.inv),
                               ny*ny, num.dist) %*% dist.weights.last, ny, ny)
    diag(dist.inv) <- 0
    
    sparseZ = Z
    if(sparse){
        ind <- which(Z==1, arr.ind=TRUE)
        sparseZ = sparseMatrix(ind[,1], j=ind[,2], x= rep(1, nrow(ind)),  dims=dim(Z))
    }
    ind <- seq.int(1, prod(dim(sparseZ)),  nrow(sparseZ))

    ## pdist = dist%*%sparseZ
    pdist = dist.inv %*% sparseZ
    if(sparse)  pdist = as(pdist, 'matrix')

    pdist0 = apply(pdist, 1, function(r) which(r==0))
    if(length(pdist0) && length(pdist0)!=ny) stop('pdist0 error')

    pdist00 = which(pdist == 0)
    
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
            ## dist.list = transform_dist(dist.original, t.max, peta[,s])
            ## dist = matrix(matrix(unlist(dist.list), ny*ny, num.dist) %*% dist.weights, ny, ny)
            ## pdist = dist %*%sparseZ
            ## pdist.list = lapply(dist.list, function(r) if(sparse) as(r %*% sparseZ, 'matrix') else r %*% sparseZ)
            
            dist.list.inv = transform_dist.inv(dist.original, t.max, peta[,s])
            dist.inv = 1/matrix(matrix(unlist(dist.list.inv),
                                       ny*ny, num.dist) %*% dist.weights.last, ny, ny)
            diag(dist.inv) <-0
            pdist = dist.inv %*%sparseZ

            pdist.list = lapply(dist.list.inv,
                                function(r)
                                    if(sparse) as(r %*% sparseZ, 'matrix')
                                    else r %*% sparseZ)
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
                    
                }
                
                ## Updating similarity matix parameter
                ywU =  if(distOnly) U0[i,]  else  y0.new[i]*(w0.new*U0[i,])
                new.eta = rEta.EB.local.multi(peta.last,
                                              pdist[i,],
                                              if(length(pdist0)) pdist0[[i]] else NULL,
                                              i,
                                              sparseZ,
                                              Z,
                                              ywU,
                                              eta_sd,
                                              sparse,
                                              dist.original,
                                              t.max,
                                              num.dist,
                                              dist.weights.last,
                                              dist.list.inv)
                
                peta.new = new.eta$eta
                peta.count = peta.count + 1*(abs(peta.new - peta.last) > tol.err)
                peta.sum = peta.sum + peta.new
                peta.last = peta.new

                ## if(new.eta$change){
                ##     pdist[i,] = new.eta$dist # very slow when using sparse assignment
                ## }
                #d.old = sapply(dist.list.inv, function(r) r[i,])
                d.old = new.eta$dist
                dist.w.new = rDist.weights(dist.weights.last,
                                           d.old,
                                           if(length(pdist0)) pdist0[[i]] else NULL,
                                           i,
                                           sparseZ,
                                           Z,
                                           ywU,
                                           weights_sd,
                                           num.dist,
                                           weights_hyper,
                                           sparse
                                           )
                dist.weights.new = dist.w.new$weights
                dist.weights.count = dist.weights.count + 1*(abs(dist.weights.new[1] - dist.weights.last[1]) >tol.err)
                dist.weights.sum = dist.weights.sum + dist.weights.new
                dist.weights.last = dist.weights.new
                
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
            weights_sd = weights_sd * exp(beta*(dist.weights.count/subItra-0.44)/log(1 +s))
            
            if(!distOnly){
                w0[,s+1] <- w0.sum/subItra
                w0.new = w0.count = w0.sum=0
                w0.last = w0[,s+1]
                y0[,s+1]<- y0.sum/subItra
                y0.new = y0.count = y0.sum=0
                y0.last = y0[,s+1]
            }
            
            peta[,s+1]<- peta.sum/subItra
            peta.new = peta.count = peta.sum =0
            peta.last = peta[,s+1]

            dist.weights[,s+1] <- dist.weights.sum/subItra
            dist.weights.new = dist.weights.count = dist.weights.sum =0
            dist.weights.last = dist.weights[,s+1]
            
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
    ## plot(1:5001, peta[2,], type = 'l')
    ## mean(peta[2,])
    ## sd(peta[2,])

    peta = peta[,-burn.in]
    dist.weights = dist.weights[, -burn.in]
    if(uncertainty)  g0 = g0[-burn.in] else g0 = NULL

    list(w = w0, y = y0, eta = peta, g = g0,
         burn.in = max(burn.in)-1,
         sd = list(w=if(distOnly) NULL else w_sd,
                   y = if(distOnly) NULL else y_sd,
                   eta= eta_sd,
                   weights = weights_sd),
         dist.weights = dist.weights)
}
