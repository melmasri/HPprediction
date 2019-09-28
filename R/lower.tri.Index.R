lower.tri.Index <-
function(n){
    ## returning the vectorized indices of the lower triangular
    ## elements of the matrix, excluding diagonal 
    unlist(sapply(1:(n-1), function(i) n*(i-1) + (i+1):n))
}
