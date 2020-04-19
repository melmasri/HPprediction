rEta.cophenetic <-
function(eta.old,tree,tree.ht,pdist.old, pd0,i, Z, ywU, eta_sd =0.01){
    eta.prop = eta.old + rnorm(length(eta.old), 0, sd = eta_sd)
    dist = 1/cophenetic(eb.phylo(tree, tree.ht, eta.prop))
    diag(dist)<-0
    pdist.new = c(crossprod(dist[i,],Z))
    no0 = which(pd0)
    if(length(no0)){
        if(length(no0)==sum(Z[i,])) likeli = -Inf else 
        likeli = sum(((log(pdist.new)- log(pdist.old))*Z[i,])[-no0])-
            sum((ywU*(pdist.new - pdist.old))[-no0]) 
    }else{
        likeli = sum((log(pdist.new)- log(pdist.old))*Z[i,] )-
            sum(ywU*(pdist.new - pdist.old))
    }
    ratio = min(1, exp(likeli));ratio
    u = (runif(1)<=ratio)
    if(u) { eta.old  = eta.prop; pdist.old = pdist.new}
    
    list (eta=eta.old, dist=pdist.old)
}
