cophFast <-
function(tree, lowerIndex, upperIndex, n, dummyMat){
    ## a fast version of phylo tree using phangorn:::coph
    ## phangorn:::coph returns a dist object.
    ## dummyMat = matrix(0, n,n)
    dummyMat[lowerIndex]<- dummyMat[upperIndex]<-(1/phangorn:::coph(tree))
    dummyMat
}
