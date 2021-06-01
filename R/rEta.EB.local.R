
rEta <-function(eta.old,
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
                kernel_func,
                kernel_name,
                ...){
    ## a faster version of rEta using faster cophenetic and sparse matrices
    change = FALSE
    if(length(no0) && length(no0)==sum(Z[i,])){
        likeli = -Inf
    }else{
        eta.prop = eta.old + eta_sd*rnorm(1)
        
        if(kernel_name =='radial' && eta.prop<=0) eta.prop = -eta.prop
        dist0 = 1/kernel_func(dist=dist.org[i,],
                              tmax=tmax,
                              eta=eta.prop,
                              param=...$param)
        dist0[i]<-0
        pdist.new = dist0%*%sZ
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

## ## ## ## ## ## ### Testing parametes
## i=2
## eta.old = peta.last
## pdist.old = pdist[i,]
## no0  = if(length(pdist0)) pdist0[[i]] else NULL
## i = i
## sZ = sparseZ
## Z = Z
## ywU = if(distOnly) U0[i,]else  y0.new[i]*(w0.new*U0[i,])
## dist.org = dist.original
## weights =                                      dist.weights.last
## dist.list = dist.list.inv
## peta.last
##
## eta_sd = c(0.1,0.1)
rEta.multi <-function(eta.old,
                      pdist.old,
                      no0,
                      i,
                      sZ,
                      Z,
                      ywU,
                      eta_sd =0.01,
                      sparse,
                      dist.org,
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
    switch_sign = sapply(dist.org, function(r) 1 *(r$kernel == 'radial'))
    if(min(switch_sign * eta.prop) < 0){
        eta.prop = eta.prop *(1 - switch_sign) - eta.prop * switch_sign
    }
    d = sapply(1:num.dist, function(j){
        kernel_func = dist.org[[j]]$kernel_func
        dd = dist.org[[j]]$dist[i,]
        tmax = dist.org[[j]]$t.max
        kernel_name = dist.org[[j]]$kernel
        param = dist.org[[j]]$param
        a = kernel_func(dist=dd,
                        eta=eta.prop[j],
                        tmax=tmax,
                        param=param
                        )
        a
    })
    d.old = sapply(dist.list, function(r)  r[i,])
    res = list()
    ord  = sample.int(num.dist,num.dist)
    k = 1
    a <- d.old
    pdist.last = pdist.old
    for (j in ord){
       ##  a <- d.old
        a[,j] <- d[,j]
        is_diff =sum(abs(d[,j] - d.old[,j])[-no0], na.rm = TRUE) > 1e-8
        pd = c(1/(a%*%weights))
        pd[i] <- 0
        pdist.new = if(sparse) (pd %*% sZ)@x else pd %*%sZ
        if(length(no0)){
            likeli = sum(((log(pdist.new)- log(pdist.last))*Z[i,])[-no0])-
                sum((ywU*(pdist.new - pdist.last))[-no0])
        }else{
            likeli = sum((log(pdist.new)- log(pdist.last))*Z[i,])-
                sum(ywU*(pdist.new - pdist.last))
        }
        if(is_diff && !is.na(likeli) && runif(1)<= min(1, exp(likeli))){
            res[[k]] <- list(eta  = eta.prop[j],
                             dist = d[,j],
                             change=TRUE)                             
            pdist.last = pdist.new
        }else{
            a[,j] <-d.old[,j]
            res[[k]]<-list(eta=eta.old[j],
                           dist=d.old[,j],
                           change=FALSE)
        }
        k = k + 1
    }
    
    j = order(ord)
    eta.new = sapply(res, function(r) r$eta)[j]
    change = any(sapply(res, function(r) r$change))
    dist.new = sapply(res, function(r) r$dist)[,j]
    pdist.new = NULL 
    list(eta=eta.new,
         pdist=c(pdist.new),
         change=change,
         dist = dist.new)

}
