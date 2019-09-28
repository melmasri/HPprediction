rEta.copheneticSuperFast <-
function(eta.old,pdist.old, no0,i, sZ, Z, ywU,
                              eta_sd =0.01, sparse, dist.org, tmax){
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
