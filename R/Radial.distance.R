Radial.distance <-function(dist, eta, ...)
{
    mc <- list(...)
    alpha = if(is.null(mc$param$alpha)) 1 else mc$param$alpha
    a = alpha*exp(-0.5*(dist/eta)^2)
    ##dist * 0
    #A = (a < tol.err) * 1
    #A * tol.err + (1-A) * a
    1/a
}
