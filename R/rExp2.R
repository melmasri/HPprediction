rExp2 <-
function(l, g, Z, Z0, Z00){
    ## Sampling from a zero-inflated Gumbel or a one-inflated Exponential
    ## Method one: P(z=0|g) = 1\delta_{s=0} + g\delta_{s>0}; dirac delta
    p = 1- exp(-l)
    unif = matrix(runif(length(l)), dim(p))
    U = 1 + 0*p
    U[-Z0] = -log(1- unif[-Z0]*p[-Z0])/(l[-Z0] + tol.err)
    aa = Z00 & (unif < g*p/(tol.err + g*p + 1-p))
    U[aa] =  -log(1 - unif[aa]*(g*p[aa] + 1-p[aa])/g)/l[aa]
    U
}
