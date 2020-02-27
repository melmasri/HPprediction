plot_Z <-
function(Z, xlab, ylab, ...){
    ## ploting interaction matrix as a binary image
    if(missing(ylab)) ylab = 'hosts'
	if(missing(xlab)) xlab = 'parasites'
    par(mar = c(5,5,1,1)+0.1)
	image(1:ncol(Z), 1:nrow(Z), t(Z[nrow(Z):1,]),
		col = c('white', 'black'), ylab=ylab, xlab=xlab,
          useRaster=TRUE,srt=45, axes=FALSE,cex.lab=3)
           

    a = 100*max(1,round(ncol(Z)*0.25 /100,0))
	axis(1, at = a*0:(ceiling(ncol(Z)/a)), cex.axis= 2)
    b = 100*max(1,round(nrow(Z)*0.25 /100,0))
    axis(2, at = b*0:(ceiling(nrow(Z)/b)),cex.axis = 2)
}
