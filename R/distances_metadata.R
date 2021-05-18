distances_metadata<-function(distances)
{
    if(!'kernel' %in% names(distances)){
        ## list of lists (mutli distance)
        A = lapply(distances,
                   function(r)
                       list(dist = if(is.phylo(r$dist)) unname(cophenetic(rescale(r$dist, 'EB', 0))/2) else unname(r$dist), 
                            kernel_func  = get_kernel(r$kernel),
                            kernel = r$kernel,
                            t.max = get.max.depth(r$dist)))
    }else{
        A = list(dist = if(is.phylo(distances$dist))
                            unname(cophenetic(rescale(distances$dist, 'EB', 0))/2) else unname(distances$dist), 
                 kernel_func  = get_kernel(distances$kernel),
                 kernel = distances$kernel,
                 t.max = get.max.depth(distances$dist))
        
    }
    return(A)
}




get_kernel <-function(kernel = 'EB'){
    if(kernel == 'EB') return(EB.distance)
    if(kernel == 'radial') return(Radial.distance)
    stop(paste0('Kernal ', kernel, ' is not supported, only EB and radial are!'))
}
