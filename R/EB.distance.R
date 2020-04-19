EB.distance <-
function(dist, tmax, a){
    if(a==0) return(2*dist)
    (2/a)*exp(a*tmax)*(1-exp(-a*dist))
}
