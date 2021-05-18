EB.distance <-
function(dist, tmax, eta, ...){
    if(eta==0){
        a = 2*dist
    }else{ 
        a = (2/eta)*exp(eta*tmax)*(1-exp(-eta*dist))
    }
    1/a
    #0 * a
}
