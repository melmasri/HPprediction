rg <-
function(Z,l){
    ## Sampling the uncertainty parameter
    ## Method one: P(z=0|g) = 1\delta_{s=0} + g\delta_{s>0}; dirac delta
    ZZ = 1*(l<1)                  # 1-inflated exponential
    ## Method 1
    M = sum(ZZ*Z)              # l<1(OR S>0) and Z=1,  N++
    N = sum((1-Z)*ZZ)          # l<1(OR S>0) and Z=0,  N-+
    rbeta(1 , N + 1, M + 1)         
}
