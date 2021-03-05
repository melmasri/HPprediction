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
        ## w.prop =0.1
        ## w.prop = weights.old[1]*exp(epsilon)
        ## w.prop = max(min(abs(w.prop),1-1e-4), 1e-4)
        ## w.prop = c(w.prop, 1-w.prop)
        ## d.old = sapply(dist.list, function(r) r[i,])
        pd.old = c(1/dist.list %*% weights.old)
        pd.old[i] <-0
        pdist.old = if(sparse) (pd.old %*% sZ)@x else pd.old %*% sZ
        
        pd = c(1/(dist.list %*% w.prop))
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

## weights.hyper = 10*rep(1,num.dist)/ num.dist
## x = seq(0, 1, 0.1)
## prop = sapply(x, function(w){
##     w.prop = c(w , 1-w)
##     ## d.old = sapply(dist.list, function(r) r[i,])
##     d.old = dist.list
##     pd.old = c(1/d.old %*% weights.old)
##     pd.old[i] <-0
##     pdist.old = if(sparse) (pd.old %*% sZ)@x else pd.old %*% sZ
##     pd = c(1/(d.old %*% w.prop))
##     pd[i] <-0
##     pdist.new = if(sparse) (pd %*% sZ)@x else pd %*%sZ
##     likeli = sum(((log(pdist.new)- log(pdist.old))*Z[i,])[-no0])-
##                 sum((ywU*(pdist.new - pdist.old))[-no0]) + c(weights.hyper %*% (log(w.prop) - log(weights.old)))
##             exp(likeli)
## })
## plot(x, prop, type = 'l')


## weights.hyper
