sample_parameter <-
function(param, MODEL,Z, tree, size = 1000, weights=NULL){
    sample_mcmc<-function(mcmc_sample,nObs,  size =1000){
        if(is.matrix(mcmc_sample))
            mcmc_sample[, sample.int(nObs, size, replace = TRUE)] else
        mcmc_sample[sample.int(nObs, size, replace = TRUE)]
    }
    if(!is.null(weights) & length(weights)!=size)
        stop('weights are not the same sampling size.')
    t.max = get.max.depth(tree)
    dist.original = unname(cophenetic(rescale(tree, 'EB', 0)))/2
      if(grepl('dist', MODEL)) {
        Y = W = 1
    }else{
        Y = sample_mcmc(param$y, ncol(param$y), size)
        W = sample_mcmc(param$w, ncol(param$w), size)
    }
    if(grepl('(dist|full)', MODEL)){
        Eta = sample_mcmc(param$eta, length(param$eta), size)
    }else Eta = 1

    
    zeroZ = which(Z>0)
    P <- 0
    for(s in 1:size){
        ## aux = sapply(1:size, function(s){
        ## setting affinity to 1 in distance model
        if(grepl('(aff|full)', MODEL)){
            YW = outer(Y[,s], W[,s])
        } else YW = 1
        ## Full or distance model
        ## Creating distance
        if(grepl('(full|dist)', MODEL)){
            distance = 1/EB.distance(dist.original, t.max, Eta[s])
            diag(distance)<-0
            distance = distance %*% Z
            distance[distance==0] <- if(grepl('dist', MODEL)) Inf else 1
        } else distance = 1
        ## models
        ## P = 1-exp(-outer(Y, W))                 # affinity model
        ## P = 1-exp(-distance)                    # distance model
        ## P = 1-exp(-YW*distance)                 # full model
        ## All in one probability matrix
        if(!is.null(weights)){
            a =  1-  exp(-YW*distance)
            Pg = a * weights[s] /(1-a + weights[s] * a)
            Pg[zeroZ] <- a[zeroZ]
            P = P + Pg
        }else{
            P = P + 1-  exp(-YW*distance)
        }
    }
    ## })
    matrix(P/size, nrow = nrow(com), ncol = ncol(com))

}
