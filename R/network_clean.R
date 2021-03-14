#' A function to clean the network and tree:
#' 
#' - converts Z to binary
#' - removes empty columns from Z
#' - confirms tree is a phylo object when using full or distance models
#' - removes tree tips that do not exist in the rows of Z
#' - removes the rows of Z with no corresponding tips in the tree
#' - orders the rows of Z according to the cophenetic function
#' 
#'
#' @param Z bipartite interaction matrix. Rows should correspond to species in the distance/tree. There should be no empty columns.
#' @param tree 'phylo' object, or an symmetric non-negative matrix of pairwise distance between the rows of Z
#' @param model.type Indicates the model to use: one of 'full', 'distance', 'affinity'
#' @param normalize logical whether to normalize the \code{tree} by maximum distance. (Default is TRUE).
#' @return Returns a list of the cleaned Z and distance objects
#'
#' @examples
#' 
#' # Simluate a tree and Z matrix
#' \dontrun{
#' tree <- rcoal(50)
#' Z <- matrix(rbinom(50*200, 1, 0.01), nrow=50, ncol=200)
#' Z <- Z[,colSums(Z)>0]
#' rownames(Z) <- tree$tip.label
#' 
#' # Clean the network and tree
#' out <- network_clean(Z, tree, model='full')
#' str(out)
#' }
#' @export
#' 
network_clean_single <-
function(Z, tree = NULL, model.type = c('full', 'distance', 'affinity'), normalize=TRUE){    
    
    require(geiger)
    require(phangorn)
    require(Matrix)

    model.type= tolower(model.type[1])

    ## General warnings are checks
    if(missing(Z)) stop('Interaction matrix is missing!')
    if(!all(range(Z)==c(0,1))) {
        warning('Z is converted to binary!', immediate. = TRUE, call.= FALSE)
        Z = 1*(Z>0)
    }
    if(grepl('aff', model.type)){
        if(!is.null(tree))
            warning('affinity model is chosen; tree can be used as a constant if parameter constant_distance = TRUE!, otherwise ignored',
                    immediate. = TRUE, call. = FALSE)
    }
    if(grepl('(dist|full)', model.type)){
        if(is.null(tree))
            stop('distance model is chosen, but tree is null!')
    }
    ## testing tree is a phylo object
    if(is.matrix(tree))
        if(!(isSymmetric(tree) & min(tree) >=0)){
            stop('tree must be a phylo object! or a symmetric matrix of non-negative pairwise distances between the rows of Z')
        }
   
    ##Phylo tree
    if(is.phylo(tree)){
        ## testing that all tips exist in Z
        if(!all(tree$tip.label %in% rownames(Z))){
            warning('not all species in tree exist in Z; missing are removed from tree!',
                    immediate.= TRUE, call. = FALSE)
            tree = drop.tip(tree, tree$tip.label[!(tree$tip.label %in% rownames(Z))])
        }
        ## Testing all names in com exist in dist
        if(!all(rownames(Z) %in% tree$tip.label)){
            warning('not all rows in Z exist tree; missing are removed from Z!',
                    immediate.= TRUE, call. = FALSE)
            Z = Z[rownames(Z) %in% tree$tip.label,]
        }
        ## removing any column of Z that has 0 interactions -- must exist here.
        if(any(colSums(Z)==0)){
            print('Z has empty columns - these have been removed!')
            Z = Z[,which(colSums(Z)>0)]
        }
        ## Testing that names in Z exist only once
        if(!all(sapply(sapply(rownames(Z),
                              function(r) which(r==tree$tip.label)), length)==1))
            stop('some row-species in Z exist more than once tree!')
        ## Making sure the order of hosts in Z and tree are the same
        aux  = cophenetic(tree)
        if(!all(rownames(aux)==rownames(Z))){
            row.order <- sapply(rownames(aux),function(r) which(r==rownames(Z)))
            print('ordering the rows of Z to match tree...')
            Z = Z[row.order,]
        }
        if(normalize){
            print('normalizing tree edges by the maximum pairwise distance!')
            tree$edge.length = tree$edge.length/(max(aux)/2)
        }
        tree = cophenetic(tree)/2
    }else{                              ## distance tree
        ## testing that all tips exist in Z
        if(!all(rownames(tree) %in% rownames(Z))){
            warning('not all row names in tree exist in Z; missing are removed from tree!',
                    immediate.= TRUE, call. = FALSE)
            aux= (rownames(tree) %in% rownames(Z))
            tree = tree[aux, aux ]
        }
        ## Testing all names in com exist in dist
        if(!all(rownames(Z) %in% rownames(tree))){
            warning('not all rows in Z exist tree; missing are removed from Z!',
                    immediate.= TRUE, call. = FALSE)
            aux = which(rownames(Z) %in% rownames(tree))
            if(length(aux)>0)
                Z = Z[aux,]
            else
                stop('No rows of Z exist in rows of tree!')
        }
        ## removing any column of Z that has 0 interactions -- must exist here.
        if(any(colSums(Z)==0)){
            print('Z has empty columns - these have been removed!')
            Z = Z[,which(colSums(Z)>0)]
        }
        ## Testing that names in Z exist only once
        if(!all(sapply(sapply(rownames(Z),
                              function(r) which(r==rownames(tree))), length)==1))
            stop('some row-species in Z exist more than once in rownames of distance!')
        
        ## Making sure the order of hosts in Z and tree are the same
        aux  = tree
        if(!all(rownames(aux)==rownames(Z))){
            row.order <- sapply(rownames(aux),function(r) which(r==rownames(Z)))
            print('ordering the rows of Z to match tree...')
            Z = Z[row.order,]
        }
        if(normalize){
            print('normalizing distance edges by the maximum pairwise distance!')
            tree = tree/max(aux)
        }
        if(min(tree) < 0){
            stop('minimum distance is negative!, input must be non-negative')
        }
    }
    list(Z = Z, distances = tree)
}

#' A function to clean the network given multiple pairwise distance matrices or phylo tree
#' 
#' - excutes \code{network_clean_single} on Z and every passed pairwise distance
#' - finds intersection species between all distances/phyo trees and keeps only the intersection
#' - orders the rows of Z according to the cophenetic function
#' 
#'
#' @param Z bipartite interaction matrix. Rows should correspond to species in the tree. There should be no empty columns.
#' @param distances a single object or a list of phylo objects and symmetric non-negative matrix of pairwise distance between the rows of Z
#' @param model.type Indicates the model to use: one of 'full', 'distance', 'affinity'
#' @param normalize  logical whether to normalize the \code{distances} by maximum distance. (Default and recommended is TRUE).
#' @return Returns a list of the cleaned Z and and list of ordered distances
#'
#' @export
#' 
network_clean <-function(Z, distances, model.type = c('full', 'distance', 'affinity'), normalize=TRUE){    

    require(geiger)
    require(phangorn)
    require(Matrix)
    
    model.type= tolower(model.type[1])
    ## General warnings are checks
    if(is.list(distances) & !is.phylo(distances)){
        objs = lapply(distances, function(r) network_clean_single(Z, r, model.type, normalize))
    }else{
        obj = network_clean_single(Z, distances, model.type, normalize)
        return (obj)
    }
    
    extract_names<-function(x) if(is.phylo(x)) x$tip.label else rownames(x)

    all.names = lapply(objs, function(r) extract_names(r$distances))
    intersect.names = names(which(table(unlist(all.names)) == length(all.names)))

    remove_tips<-function(dist, names){
        if(is.phylo(dist)){
            ## testing that all tips exist in Z
            if(!all(dist$tip.label %in% names)){
                warning('not all species exist in all distnaces; missing are removed from tree!',
                        immediate.= TRUE, call. = FALSE)
                dist = drop.tip(dist, dist$tip.label[!(dist$tip.label %in% names)])
                if(normalize){
                    aux= cophenetic(dist)
                    print('normalizing tree edges by the maximum pairwise distance!')
                    dist$edge.length = dist$edge.length/(max(aux)/2)
                }
            }
        }else{
            if(!all(rownames(dist)  %in% names)){
                warning('not all species exist in all distances; missing are removed from distance!',
                        immediate.= TRUE, call. = FALSE)
                aux= (rownames(dist) %in% names)
                dist = dist[aux, aux ]
                
                if(normalize){
                    print('normalizing distance edges by the maximum pairwise distance!')
                    dist = dist/max(dist)
                }
            }
        }
        return(dist)
    }
    ## keeping only intersection rows/hosts and normalizing
    distances = lapply(objs, function(r) remove_tips(r$distances, intersect.names))
    
    ## keeping only species names in the intersection in Z 
    Z = objs[[1]]$Z
    if(!all(rownames(Z) %in% intersect.names)){
            warning('not all rows in Z exist in the intersection of names in distances; missing are removed from Z!',
                    immediate.= TRUE, call. = FALSE)
            aux = which(rownames(Z) %in% intersect.names)
            if(length(aux)>0)
                Z = Z[aux,]
    }
    ## removing any column of Z that has 0 interactions -- must exist here.
    if(any(colSums(Z)==0)){
        print('Z has empty columns - these have been removed!')
        Z = Z[,which(colSums(Z)>0)]
    }
    
    ## ordering all quantities based on the phylo tree if possible,
    ## Ordering all distances and rows of Z
    order_dist <-function(d, names1, twosides = TRUE) {
        ## we order everything based on the phylo distnace
        #+ TODO: how to order the phylo as well based on list of names
        if(!is.phylo(d)){
            aux = pmatch(names1, rownames(d))
            if(twosides)
                return (d[aux,aux])
            else return(d[aux,])
        }
        return(d)
    }
    ## finding a ordering
    order.names = unlist(sapply(distances, function(r) if(is.phylo(r)) rownames(cophenetic(r))))
    if(is.null(order.names)) order.names = intersect.names
    
    distances.ordered = lapply(distances, function(r) order_dist(r, order.names))
    Z   = order_dist(Z, order.names, FALSE)
    
    list(Z = Z, distances = distances.ordered)
}
