row_extract <-
function(m, i, ind){
    ## extracting row i from matrix m based on ind
    ## m: matrix of package Matrix
    ## i: is the row number
    ## ind: is the index of the first cell of each column
    ## ind can be constructed as
    ## seq.int(1, prod(dim(m)),  nrow(m))
    
    if(class(m)=='dgeMatrix')
        m@x[ind + (i-1)] else m[ind + (i-1)]
}
