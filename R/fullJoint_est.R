fullJoint_est <-
function(Z, iter = 10, uncertainty = FALSE, distance = NULL, ...){
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

    if(!is.null(el$constant_distance) & !is.null(distance)){
        if(!isSymmetric(distance) & nrow(Z)!=nrow(distance))
            stop('distance must be symmetrics with number of rows equal to the rows of Z!')
        pdist = distance %*% Z
    }else{
        pdist = 1
    }
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
                    save(y0,w0,g0, file='snapshot.RData')
            }
            ## Updatting latent scores
            U0 <- rExp(outer(y0[,s],w0[, s]) * pdist)
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
    list(w = w0, y = y0, g = g0, burn.in = max(burn.in)-1, distance = pdist,
         sd = list(w=w_sd, y = y_sd))
}
