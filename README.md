# Repository description

The `ICM` branch fits the model using the "Iterative Conditional Modelling" approach, discussed Besag(74) Section 6.1, which is called the "Coding" method.

## Quick intro

The main function is `network_est()` in file `networkMCMC.R`

`network_est` returns a a list of input objects, the binary interaction matrix and the phylogeny tree in some cases, and the MCMC samples of the specified parameters interest.

## Usage

`network_est(Z, slices = 10, tree = NULL, model.type = c('full', 'distance', 'affinity'), uncertainty = FALSE, ... )`

## Arguments
    + `Z`:  an H x J binary or count matrix of interactions between two sets of species, H and J. If `Z` is count it will be converted to binary.
    + `slices`:" the number of slices or MCMC samples to run. In the case of the full model or the phylogeny-only model a single sample/itration is the average of H samples from H conditional distributions, one for each row of `Z`; as in the ICM model. For the affinity-only model sampling from the full joint is possible; see 'Details' for more information.
    + `model.type`: either distance, affinity, or full to combine both distance and affinity. Default is the full.
    + `uncertainty`: whether to sample an uncertainty parameter of not. See 'Details' for more information.
    + `...`: optional arguments to the lower layer MCMC sampling algorithm:
        + `y`: a single number of a vector of size H specifying the initial value for row affinities, default is 1. 
        + `a_y`: alpha hyperparameter of the Gamma priors of row affinities, default is 1. The beta parameter is set to 1. 
        + `y_sd`: initial standard deviation for the proposal distribution (Gaussian) for the row parameters, default is 0.2. 
        + `w`: a single number of a vector of size J specifying the initial value for column affinities, default is 1. 
        + `a_w`: alpha hyperparameter of the Gamma priors of column affinities, default is 1. The beta parameter is set to 1.
        + `w_sd`: initial standard deviation for the proposal distribution (Gaussian) for the column parameters, default is 0.2.
        + `eta`: initial value for tree transformation parameter under the `EB` model, default is 0.
        + `eta_sd`: initial standard deviation for the proposal distribution (Gaussian) for the column parameters, default is 0.005.
        + `burn.in`:  proportion of the initial samples to discard (0 to 1), default is 0.5.
        + `batch.size`: update frequency of the Adaptive MH for `_sd` parameters, default every 50 iterations. This applies to affinity-only model.
        + `beta`:  tuning parameter for updating the Adaptive MH (>=0, 0 for no tuning, larger means more aggressive tuning).


## Details

The function `network_est()` returns a list that includes a list names `param` of posterior samples for each estimated parameter with burn-in removed. In addition the input data; `Z` and `tree` for `full` and `distnace` models, and `Z` only for `affinity` model. Note that `network_est()` does initial cleaning of `Z` and `tree` to conform to the needed input. Such as:
    + converting `Z` to binary;
    + removing row-species in `Z` that are not in `tree`;
    + removing empty rows and columns from `Z`;
    + removing tips in `tree` that do not correspond to row-species in `Z`;
    + left ordering of `Z`.
    
For convenience, row parameters (gammas), and column parameters (rhos), are indicated by the variables `y` and `w`, respectively, in the output. If using the `distance` model or `full`, this will also include posterior samples for the phylogenetic transformation parameter `eta`. The uncertainty parameter is denoted by `g`, and returned only when `uncertainty=TRUE`, otherwise `NULL`.

Specifying inital values for affinity parameters and related options in `...` is only relavent for the `full` or `affinity` model, otherwise they are ignored. Intial values for the `eta` parameter only apply in the `full` or `distance` model.
    
`network_est()` runs two types of upper layer MCMC sampling functions. It runs and MCMC using the full joint distribution when the `affinity` model is used, otherwise it runs an iterated conditional modes (ICM) MCMC sampler when the `full` or `distance` models are used. A mode of the ICM is the conditional joint distribution of a row on all other rows of `Z`. 
    
## Examples



    
