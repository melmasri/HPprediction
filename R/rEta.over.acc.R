rEta.over.acc <-
function(eta.old, pdist.old, pdist0,pdist00, sZ, Z, ywU,
                        eta_sd =0.01, sparse, dist.org, tmax){
    ## a faster version of rEta using faster cophenetic and sparse matrices
    eta.prop = eta.old + eta_sd*rnorm(1)
    dist = 1/EB.distance(dist.org, tmax, eta.prop)
    diag(dist)<-0
    pdist.new =  dist%*%sZ
    if(sparse) pdist.new = as(pdist.new, 'matrix') 
    pdist.new[pdist00]<-1
    acc = sapply(1:nrow(Z), function(i){
        no0  = if(length(pdist0)) pdist0[[i]] else NULL
        if(length(no0)){
            if(length(no0)==sum(Z[i,])) likeli = -Inf else 
            likeli = sum(((log(pdist.new[i,])- log(pdist.old[i,]))*Z[i,])[-no0])-
                sum((ywU[i,]*(pdist.new[i,] - pdist.old[i,]))[-no0]) 
        }else{
            likeli = sum((log(pdist.new[i,])- log(pdist.old[i,]))*Z[i,] )-
                sum(ywU[i,]*(pdist.new[i,] - pdist.old[i,]))
        }
        likeli
    })
    ##acc
    acc = mean(pmin(1, exp(acc)))
    if(!is.nan(acc) && runif(1)<= acc)
        return(list(eta = eta.prop,dist= pdist.new,change=TRUE))
    else return(list(eta = eta.old,change=FALSE))
    
}
