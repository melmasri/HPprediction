## i=2
## weights.old=dist.weights.last
## dist.list = d.old
## no0=if(length(pdist0)) pdist0[[i]] else NULL
## sZ=sparseZ
## Z=Z
## weights_sd=weights_sd
## weights.hyper=weights_hyper

rDist.weights <-function(weights.old,
                         dist.list,
                         no0,
                         i,
                         sZ,
                         Z,
                         ywU,
                         weights_sd,
                         num.dist,
                         weights.hyper,
                         sparse
                         ){
    change = FALSE
    if(length(no0) && length(no0)==sum(Z[i,])){
        likeli = -Inf
    }else{
        epsilon = weights_sd * rnorm(num.dist)
        w.prop = weights.old * exp(epsilon)/sum(weights.old * exp(epsilon))

        pd.old = c(dist.list %*% weights.old)
        pd.old[i] <-0
        pdist.old = if(sparse) (pd.old %*% sZ)@x else pd.old %*% sZ
        pd = c(dist.list %*% w.prop)
        pd[i] <-0
        pdist.new = if(sparse) (pd %*% sZ)@x else pd %*%sZ
        
        if(length(no0)){
            likeli = sum(((log(pdist.new)- log(pdist.old))*Z[i,])[-no0])-
                sum((ywU*(pdist.new - pdist.old))[-no0]) + c(weights.hyper %*% (log(w.prop) - log(weights.old)))
        }else{
            likeli = sum((log(pdist.new)- log(pdist.old))*Z[i,] )-
                sum(ywU*(pdist.new - pdist.old)) + c(weights.hyper %*% (log(w.prop) - log(weights.old)))
        }
    }
    if(!is.na(likeli) && runif(1)<= min(1, exp(likeli)))
    {weights.old  = w.prop;change=TRUE}
    
    list (weights =weights.old, change=change)
}
