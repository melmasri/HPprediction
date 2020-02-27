upper.tri.Index <-
function(n){
    ## returning the vectorized indices of the upper triangular
    ## elements of the matrix, excluding diagonal 
    b = matrix(0, n ,n )
    lowerIndex = lower.tri.Index(n)
    b[lowerIndex]<-lowerIndex
    b =t(b)
    a = which(upper.tri(b))
    a[order(b[upper.tri(b)])]
}
