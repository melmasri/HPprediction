# Repository description

The `ICM` branch fits the model using the "Iterative Conditional Modelling" approach, discussed Besag(74) Section 6.1, which is called the "Coding" method.

# Usage

The main function is `network_est()` in file `networkMCMC.R`


```
network_est returns a a list of input objects, the binary interaction matrix and the phylogeny tree in some cases, and the MCMC samples of the specified parameters interest.

Usage

network_est(Z, slices = 10, tree = NULL, model.type = c('both', 'distance', 'affinity'), uncertainty = FALSE, ... )

Arguments
Z 	

an H x J binary or count matrix of interactions between two sets of species, H and J. If Z is count it will be converted to binary.


slices

the number of slices or MCMC samples to run. In the case of the full model or the phylogeny-only model a single sample/itration is the average of H samples from H conditional distributions, one for each row of Z; as in the ICM model. For the affinity-only model sampling from the full joint is possible; see 'Details' for more information.


model.type

either the full model (both), phylogeny-only (distance) or the affinity-only (affinity) models. Default is the full model (both)


uncertainty

whether to sample an uncertainty parameter of not. See 'Details' for more information.

... 	

optional arguments to the lower layer MCMC sampling algorithm. See 'Details' for more information.

Details

```



      

    
