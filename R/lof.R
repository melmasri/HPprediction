#' A function to left order a matrix based on first occurrence per row
#'
#' @param Z a binary matrix
#' @param indices a logical parameter whether to return the column order (`TRUE`) or to return the left-ordered matrix. Defautl is `TRUE`.
#' 
#' @description Given a binary matrix Z with fixed rows, `lof` left orders the matrix sequentially based on first occurrence of each row
#' 
#' 
#' @return return the left-ordered matrix, or a vector of indexing the order of columns
#' @examples
#' \donttest{
#' Z <- matrix(rbinom(10*30, 1, 0.1), nrow=10, ncol=30)
#' Zord = lof(Z)
#' }
#' 
#' @export
#' 
lof <-
function(Z, indices = FALSE){
    if(min(range(Z))<0) stop('Range is less that 0.')

    active_col <- apply(Z,1,function(r) which(as.vector(r)>0))
	bank = active_col[[1]]
	for(i in 1:nrow(Z)){
			a = setdiff(active_col[[i]], bank)
			bank=c(bank,a)
		}
    if(indices) bank else  Z[,bank]
}
