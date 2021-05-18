Radial.distance <-function(dist, eta, alpha=1, ...)
{
    tol.err = 1e-7
    a = alpha*exp(-0.5*(dist/eta)^2)
    ##dist * 0
    A = (a < tol.err) * 1
    A * tol.err + (1-A) * a
}
