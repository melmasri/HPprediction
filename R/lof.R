lof <-
function(Z, indices = FALSE){
    ## Given a binary matrix Z. Where the rows is fixed
    ## A function that left orders the matrix sequentially from row 1 to n
    ## based on first appearance of columns.
    if(min(range(Z))<0) stop('Range is less that 0.')

    active_col <- apply(Z,1,function(r) which(as.vector(r)>0))
	bank = active_col[[1]]
	for(i in 1:nrow(Z)){
			a = setdiff(active_col[[i]], bank)
			bank=c(bank,a)
		}
    if(indices) bank else  Z[,bank]
}
