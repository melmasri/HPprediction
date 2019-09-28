topPairs <-
function(P,Z,topX=20){
    ## Returning pairs with highest posterior probability
    require(reshape2)
    P[Z>0]<--1
    aux =   melt(P)
    aux = aux[order(aux$value,decreasing=TRUE),]
    colnames(aux)<-c('Hosts', 'Parasites', 'p')
    aux[1:topX,]
}
