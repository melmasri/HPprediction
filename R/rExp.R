rExp <-
function(l,a=1){
    ## Sampling from a truncated exponential distribution
    tol.err = 1e-4
    unif = runif(length(l))
	-log(1- unif*(1-exp(-l*a)))/(l + tol.err)
}
