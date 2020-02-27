network_est <-
    function(Z, slices = 10, tree = NULL, model.type = c('full', 'distance', 'affinity'),
             uncertainty = FALSE,
             sparse=TRUE, ... ){
        require(geiger)
        require(phangorn)
        require(Matrix)
        require(methods)
        ## Running options:
        ## 1 - Full to combined 2 and 3;
        ## 2 - Affinity-only model (affinity);
        ## 3 - distance-only model (distance).
        model.type= tolower(model.type[1])
        ## General warnings are checks
        if(missing(Z)) stop('Interaction matrix is missing!')
        if(slices==0)  stop('no. of slices cannot be 0!')

        cleaned = network_clean(Z, tree, model.type)
        Z = cleaned$Z
        tree = cleaned$tree

        ## For affinity-only model 
        if(grepl('aff', model.type)){
            ## sparse option is not used for affinity
            print('Running affinity model...')
            param = fullJoint_est(Z, iter = slices, uncertainty = uncertainty, ...)
        }
        
        ##  Full and distance model
        el <-list(...)
        if(is.null(el$experimental.ICM)){
            if(grepl('(dist|full)', model.type)){
                ## Running the MCMC
                print(paste0('Running ',
                             ifelse(grepl('dist', model.type),
                                    'distance model...', 'full model...')))
                param  = ICM_est(unname(Z),tree,slices, distOnly = grepl('dist', model.type),
                                 uncertainty = uncertainty, sparse=sparse, ...)
            }
        }else{
            if(grepl('(dist|full)', model.type)){
                ## Running the MCMC
                print(paste0('Running ',
                             ifelse(grepl('dist', model.type),
                                    'distance model...', 'full model...')))
                warning('Running an experimental ICM procedure!', immediate. = TRUE, call.= FALSE)
                param  = ICM_est_over_acc(unname(Z),tree,slices, distOnly = grepl('dist', model.type),
                                          uncertainty = uncertainty, sparse=sparse, ...)
            }
            
        }
        list(param =param , Z = Z, tree=tree)
    }
