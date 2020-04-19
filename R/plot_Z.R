#' A function to plot a bipartite interaction matrix
#'
#' @param Z bipartite interaction matrix
#' @param xlab labels for the x axis (default is 'parasites')
#' @param ylab labels for the y axis (default is 'hosts')
#' @param ... Additional parameters that control graphing parameters
#' 
#' @description
#' 
#' This function generates a raster plot of the binary interaction matrix 'Z' using the 'image' function in the 'graphics' package.
#'
#' @examples
#' 
#' # Simluate a Z matrix and plot it
#' \dontrun{
#' Z <- matrix(rbinom(50*200, 1, 0.01), nrow=50, ncol=200)
#' Z <- Z[,colSums(Z)>0]
#' plot_Z(tree)
#' }
#' @export
#' 
plot_Z <-
function(Z, xlab, ylab, tickMarks=100, ...){
    ## ploting interaction matrix as a binary image
    if(missing(ylab)) ylab = 'hosts'
    if(missing(xlab)) xlab = 'parasites'
    # par(mar = c(5,5,1,1)+0.1)
    image(1:ncol(Z), 1:nrow(Z), t(Z[nrow(Z):1,]),
        col = c('white', 'black'), ylab=ylab, xlab=xlab,
          useRaster=TRUE,srt=45, axes=FALSE, ...)
           
    a = tickMarks*max(1,round(ncol(Z)*0.25 /100,0))
    axis(1, at = a*0:(ceiling(ncol(Z)/a)), ...)
    b = tickMarks*max(1,round(nrow(Z)*0.25 /100,0))
    axis(2, at = b*0:(ceiling(nrow(Z)/b)), ...)
}
