#' A function to list the top predicted interactions from a model
#'
#' @param P the estimated posterior probability interaction matrix
#' @param Z the observed bipartite interaction matrix
#' @param topX number of interactions to return (default is 10)
#'  
#' @description
#' 
#' This pulls out the top predicted interactions from a model
#'  
#' @return 
#' 
#' Returns a datafame of Host, Parasite, and p(interaction) with length determined by topX (default is 10)
#' 
#' @examples
#' \dontrun{
#' topPairs(P, Z, topX=15)
#' }
#' @export
#' 
topPairs <-function(P, Z, topX=10){
    P[Z>0]<--1
    aux =   melt(P)
    aux = aux[order(aux$value,decreasing=TRUE),]
    colnames(aux)<-c('Hosts', 'Parasites', 'p')
    aux[1:topX,]
}

