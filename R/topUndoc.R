#' A function to list the undocumented interactions with the highest posterior predictions
#'
#' @param P the estimated posterior probability interaction matrix
#' @param Z the observed bipartite interaction matrix
#' @param topX number of interactions to return (default is 10)
#'  
#' @description
#'  
#' This pulls out the top undocumented interactions with the highest posterior prediction.
#'  
#' @return 
#' Returns a datafame of Host, Parasite, and p(interaction) with length determined by topX (default is 10)
#'
#' @examples
#' 
#' topUndoc(P, Z, topX=15)
#'  
#' @export
#' 
topUndoc<-function(P, Z, topX=10){
    ## Plotting the top interactions that were not documented in the original data
    require(reshape2)
    P[Z==1]<--1
    rownames(P) <- rownames(Z)
    colnames(P) <- colnames(Z)
    aux = melt(P)
    aux = aux[order(aux$value,decreasing=TRUE),]
    colnames(aux)<-c('Host', 'Parasite', 'p(interaction)')
    aux
    aux[1:topX,]
}
