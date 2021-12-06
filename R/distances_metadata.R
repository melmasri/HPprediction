distances_metadata<-function(distances)
{
    if(!'kernel' %in% names(distances)){
        ## list of lists (mutli distance)
        A = lapply(distances, distances_metadata_single)
    }else{
        ## a single distance
        A = distances_metadata_single(distances)
    }
    return(A)
}

distances_metadata_single<-function(D)
{
    A = list(dist = if(is.phylo(D$dist))
                        unname(cophenetic(rescale(D$dist, 'EB', 0))/2) else unname(D$dist), 
             kernel_func  = get_kernel(D$kernel),
             kernel = D$kernel,
             t.max = get.max.depth(D$dist),
             param = D$param
             )
return(A)
    
}



get_kernel <-function(kernel = 'EB', ...){
    if(kernel == 'EB') return(EB.distance)
    if(kernel == 'radial') return(Radial.distance)
    if(kernel == 'identity') return(Identity.distance)
    stop(paste0('Kernal ', kernel, ' is not supported, only EB, radial and identity are!'))
}
