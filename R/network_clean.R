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
#' @param Z bipartite interaction matrix. Rows should correspond to species in the tree. There should be no empty columns.
#' @param tree 'phylo' object
#' @param model.type Indicate model to use: one of 'full', 'distance', 'affinity'
#'
#' @return Returns a list of the cleaned Z and tree objects
#'
#' @examples
#' 
#' # Simluate a tree and Z matrix
#' 
#' require(ape)
#' tree <- rcoal(50)
#' Z <- matrix(rbinom(50*200, 1, 0.01), nrow=50, ncol=200)
#' Z <- Z[,colSums(Z)>0]
#' rownames(Z) <- tree$tip.label
#' 
#' # Clean the network and tree
#' out <- network_clean(Z, tree, model='full')
#' str(out)
#' 
#' @export
#' 
network_clean <-
function(Z, tree = NULL, model.type = c('full', 'distance', 'affinity')){    

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
    if(any(colSums(Z)==0)){
        stop('Z has empty columns, please remove!', immediate. = TRUE, call.= FALSE)
        Z = Z[,which(colSums(Z)>0)]
    }
    if(grepl('aff', model.type)){
        if(!is.null(tree))
            warning('affinity model is chosen; ignoring tree!',
                    immediate. = TRUE, call. = FALSE)
    }
    if(grepl('(dist|full)', model.type)){
        if(is.null(tree))
            stop('distance model is chosen, but tree is null!')
        ## testing tree is a phylo object
        if(!is.phylo(tree))
            stop('tree must be a phylo object!')

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
        if(max(range(tree$edge.length))>1){
            print('normalizing tree edges by the maximum pairwise distance!')
            tree$edge.length = tree$edge.length/(max(aux)/2)
        }
    }
    list(Z = Z, tree = tree)
}
