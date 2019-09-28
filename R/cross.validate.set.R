cross.validate.set <-
function(Z, rate= 0.2){
    ##  Retunrs a Z with a percentage of one's turned 0
    ## only cells with rows of more than one interaction and columns with more than 2
    if(max(range(Z))>1) Z[Z>0]<-1
    Zo = Z
    pairs = which(Z==1, arr.ind=T)
    n = floor(rate*nrow(pairs))
    sampled = rep(0, nrow(pairs))
    repeat{
        r = sample(which(sampled==0),1)
        if(sum(Z[pairs[r,1],])>1 & sum(Z[,pairs[r,2]])>2){
            Z[pairs[r,1], pairs[r,2]]<-0
            sampled[r]<-1
            n = n -1
        }
        if(n==0) break
    }
    print(sprintf("Actual cross-validation rate is %0.3f" ,sum(Zo-Z)/sum(Zo)))
    Z
}
