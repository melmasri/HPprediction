#' A function to run an MCMC sampler based on the Iterative Conditional Modes (ICM) algorithm:
#'
#' @param Z bipartite interaction matrix
#' @param distances a 'phylo' object representing rows of Z, or a matrix of pairwise non-negative distances of rows of Z, or a list of such objects
#' @param slices The total number of samples to take
#' @param model.type Indicate model to use: one of 'full', 'distance', 'affinity'
#' @param uncertainty Indicate whether to use model variant that accounts for uncertainty in unobserved interactions.
#' @param sparse Indicate whether to use the sparseMatrix function from the Matix package. This option is advised for large matrices.
#' @param ... Additional parameters that are passed to fullJoint_est (for the 'affinity' model), and ICM_est (for the 'full' and 'distance' models) - see Description
#' 
#' @description
#' 
#' Parameter set-up
#' 
#' `y` is the initial value for the host affinity parameters
#' `w` is the initial value for the parasite affinity parameters
#' `a_y` is the shape parameter for the gamma prior on the host affinity parameters
#' `a_w` is the shape parameter for the gamma prior on the parasite affinity parameters
#' `b_y` is the rate parameter for the gamma prior on the host affinity parameters
#' `b_w` is the rate parameter for the gamma prior on the parasite affinity parameters
#' `y_sd` is the standard deviation of the proposal distribution for the host affinity parameters (normal distribution with mean zero and standard deviation y_sd).
#' `w_sd` is the standard deviation of the proposal distribution for the parasite affinity parameters (normal distribution with mean zero and standard deviation w_sd). 
#'
#' `eta` is the initial value for the tree (or distance) scaling parameter eta ('distance' and 'full' models)
#' `eta_sd` is the standard deviation of the proposal distribution for the eta parameter (normal distribution with mean zero and standard deviation eta_sd)
#' 
#' Burn-in set-up
#' 
#' `burn.in` sets the proportion of initial slices to discard from the posterior (default is 0.5)
#'
#' Metropolis Hastings adapation set-up
#'
#' `batch_size` is the Metropolis Hastings adaptiveness for updating hyperparameters. Indicated by a number of samples, this is how often the hyperparameters are updated during sampling. This is only used for the 'affinity' model, for all other models the batch size is set to 1. 
#' `beta` is the adaptation parameter for the Metropolis Hastings algorithm.
#' 
#' @return Returns a list of the posterior samples for each of the estimated parameters.
#' `y` are the host affinity parameters (for the 'affinity' and 'full  models')
#' `w` are the parasite affinity parameters (for the 'affinity' and 'full  models')
#' `eta` is the phylogeny scaling paramers (for 'distance' and 'full' models)
#' `g` is the uncertainty parameter (if uncertainty=TRUE)
#' `burn.in` is the number of slices discarded for burn-in
#' `sd` is a list of the standard deviation parameters (representing y_sd, w_sd, and eta_sd according to the model chosen)
#'  
#'
#' @examples
#' 
#' # Simluate a tree and Z matrix
#'\dontrun{ 
#' tree <- rcoal(50)
#' Z <- matrix(rbinom(50*200, 1, 0.01), nrow=50, ncol=200)
#' Z <- Z[,colSums(Z)>0]
#' rownames(Z) <- tree$tip.label
#' 
#' # Clean the network and tree
#' out <- network_est(Z, tree, model='full')
#' str(out)
#' }
#' @export
#' 
network_est <-
    function(Z, distances = NULL, slices = 10, model.type = c('full', 'distance', 'affinity'),
             uncertainty = FALSE,
             sparse=TRUE, ... ){
        require(geiger)
        require(phangorn)
        require(Matrix)
        require(methods)
        model.type= tolower(model.type[1])
        ## General warnings are checks
        if(missing(Z)) stop('Interaction matrix is missing!')
        if(slices==0)  stop('no. of slices cannot be 0!')
        
        
        cleaned = network_clean(Z, distances, model.type)
        Z = cleaned$Z
        distances = cleaned$distances
        num.distances = 1
        ## For affinity-only model 
        if(grepl('aff', model.type)){
            ## sparse option is not used for affinity
            print('Running affinity model...')
            param = fullJoint_est(unname(Z), iter = slices, uncertainty = uncertainty, distance = distances, ...)
        }
        
        ##  Full and distance model
        el <-list(...)
        if(grepl('(dist|full)', model.type)){
            if(is.null(el$experimental.ICM)){
                num.distances = if(is.list(distances) &!is.phylo(distances)) length(distances) else 1
                ## Running the MCMC
                if(num.distances  > 1){
                       print(paste0('Running multiple distances',
                                 ifelse(grepl('dist', model.type),
                                        'distance model...', 'full model...')))
                       
                    param  = ICM_est_multidistance(unname(Z),distances,slices, distOnly = grepl('dist', model.type),
                                                   uncertainty = uncertainty, sparse=sparse, ...)
                       
                }else{
                ## Running the MCMC
                    print(paste0('Running ',
                                 ifelse(grepl('dist', model.type),
                                        'distance model...', 'full model...')))
                    param  = ICM_est(unname(Z),distances,slices, distOnly = grepl('dist', model.type),
                                     uncertainty = uncertainty, sparse=sparse, ...)
                }
            }else{
                warning('Running an experimental ICM procedure!', immediate. = TRUE, call.= FALSE)
                param  = ICM_est_over_acc(unname(Z),
                                          distances, slices, distOnly = grepl('dist', model.type),
                                          uncertainty = uncertainty, sparse=sparse, ...)
                
            }
        }
        list(param=param , Z = Z, distances=distances, slices = slices, model.type = model.type, num.distances = num.distances)
    }
