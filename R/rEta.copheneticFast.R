rEta.copheneticFast <-
function(eta.old,tree,tree.ht,pdist.old, no0,i, sZ, Z, ywU,
                              eta_sd =0.01, a, b,nr, ind, sparse, dummyMat){
    ## a faster version of rEta using faster cophenetic and sparse matrices
    change = FALSE
    eta.prop = eta.old + eta_sd*rnorm(1)
    dist = cophFast(eb.phylo(tree, tree.ht, eta.prop), a, b,nr, dummyMat)
    pdist.new = dist[i,]%*%sZ
    if(sparse)
        pdist.new= pdist.new@x
    if(length(no0)){
        if(length(no0)==sum(Z[i,])) likeli = -Inf else 
        likeli = sum(((log(pdist.new)- log(pdist.old))*Z[i,])[-no0])-
            sum((ywU*(pdist.new - pdist.old))[-no0]) 
    }else{
        likeli = sum((log(pdist.new)- log(pdist.old))*Z[i,] )-
            sum(ywU*(pdist.new - pdist.old))
    }
    if(!is.na(likeli) && runif(1)<= min(1, exp(likeli)))
        { eta.old  = eta.prop; pdist.old = c(pdist.new);change=TRUE}
    
    list (eta=eta.old, dist=pdist.old, change=change)
}
