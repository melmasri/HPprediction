
rEta.EB.local <-function(eta.old,
                         pdist.old,
                         no0,
                         i,
                         sZ,
                         Z,
                         ywU,
                         eta_sd =0.01,
                         sparse,
                         dist.org,
                         tmax){
    ## a faster version of rEta using faster cophenetic and sparse matrices
    change = FALSE
    if(length(no0) && length(no0)==sum(Z[i,])){
        likeli = -Inf
    }else{
        eta.prop = eta.old + eta_sd*rnorm(1)
        dist = 1/EB.distance(dist.org[i,], tmax, eta.prop)
        dist[i]<-0
        pdist.new = dist%*%sZ
        if(sparse)
            pdist.new= pdist.new@x
        if(length(no0)){
            likeli = sum(((log(pdist.new)- log(pdist.old))*Z[i,])[-no0])-
                sum((ywU*(pdist.new - pdist.old))[-no0]) 
        }else{
            likeli = sum((log(pdist.new)- log(pdist.old))*Z[i,] )-
                sum(ywU*(pdist.new - pdist.old))
        }
    }
    
    if(!is.na(likeli) && runif(1)<= min(1, exp(likeli)))
        { eta.old  = eta.prop; pdist.old = c(pdist.new);change=TRUE}
    
    list (eta=eta.old, dist=pdist.old, change=change)
}

## ## ### Testing parametes
## i=2
## eta.old = peta.last
## pdist.old = pdist[i,]
## no0  = if(length(pdist0)) pdist0[[i]] else NULL
## i = i
## sZ = sparseZ
## Z = Z
## ywU = if(distOnly) U0[i,]else  y0.new[i]*(w0.new*U0[i,])
## dist.org = dist.original
## tmax = t.max
## weights = dist.weights
## dist.list = dist.list.inv
## ## ## ##

rEta.EB.local.multi <-function(eta.old,
                               pdist.old,
                               no0,
                               i,
                               sZ,
                               Z,
                               ywU,
                               eta_sd =0.01,
                               sparse,
                               dist.org,
                               tmax,
                               num.dist,
                               weights,
                               dist.list){
    ## a faster version of rEta using faster cophenetic and sparse matrices
    if(length(no0) && length(no0)==sum(Z[i,]))
        return(list (eta=eta.old,
                     pdist=pdist.old,
                     change=FALSE,
                     dist = sapply(dist.list, function(r){
                         r[i,]
                     }) ))
    
    eta.prop = eta.old + eta_sd*rnorm(num.dist)
    d = sapply(1:num.dist, function(j){
        EB.distance(dist.org[[j]][i,], tmax[[j]], eta.prop[j])
        ##d[i]<-0
        ## pd = d %*% sZ
        ## if(sparse)  pd@x else  pd
                                        #d
    })
    
    d.old = sapply(dist.list, function(r){
        r[i,]
    })


    ## res = lapply(1:num.dist, function(j){
    ##     a =d.old
    ##     a[,j] <- d[,j]
    ##     pdist.new = a %*% weights
    ##     if(length(no0)){
    ##         likeli = sum(((log(pdist.new)- log(pdist.old))*Z[i,])[-no0])-
    ##             sum((ywU*(pdist.new - pdist.old))[-no0]) 
    ##     }else{
    ##         likeli = sum((log(pdist.new)- log(pdist.old))*Z[i,] )-
    ##             sum(ywU*(pdist.new - pdist.old))
    ##     }
    ##     if(!is.na(likeli) && runif(1)<= min(1, exp(likeli)))
    ##         return(list(eta  = eta.prop[j], dist = d[,j],change=TRUE))
    ##     return(list(eta=eta.old[j], dist=d.old[,j], change=FALSE))
    ## })
    
    ## eta.new = sapply(res, function(r) r$eta)
    ## eta.new[2] = eta.new[1]
    ## change = any(sapply(res, function(r) r$change))
    ## pdist.new = sapply(res, function(r) r$dist) %*% weights

    a = d.old
    res = list()
    ord  = sample.int(num.dist,num.dist)
    k = 1
    for (j in ord){
        a[,j] <- d[,j]
        pd = c(1/(a%*% weights))
        pd[i] <- 0
        pdist.new = if(sparse) (pd %*% sZ)@x else pd %*%sZ
        
        if(length(no0)){
            likeli = sum(((log(pdist.new)- log(pdist.old))*Z[i,])[-no0])-
                sum((ywU*(pdist.new - pdist.old))[-no0]) 
        }else{
            likeli = sum((log(pdist.new)- log(pdist.old))*Z[i,] )-
                sum(ywU*(pdist.new - pdist.old))
        }
        if(!is.na(likeli) && runif(1)<= min(1, exp(likeli))){
            res[[k]] <- list(eta  = eta.prop[j], dist = d[,j],change=TRUE)
        }else{
            a[,j] <-d.old[,j]
            res[[k]]<-list(eta=eta.old[j], dist=d.old[,j], change=FALSE)
        }
        k = k + 1
    }

    j = order(ord)
    eta.new = sapply(res, function(r) r$eta)[j]
    change = any(sapply(res, function(r) r$change))
    dist.new = sapply(res, function(r) r$dist)[,j]
    
    ## pd = c(1/(dist.new %*%weights))
    ##pd[i] <-0
    pdist.new = NULL 
    
    list(eta=eta.new, pdist=c(pdist.new), change=change, dist = dist.new)

}
