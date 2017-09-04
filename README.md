`R` code of the implementation of the model in

Elmasri, M., Farrell, M., and Stephens, D. A. (2017). _A hierarchical Bayesian model for predicting ecological interactions using evolutionary relationships._ [arXiv](https://arxiv.org/abs/1707.08354).

## Quick intro

The main function is `network_est()` in file `network_MCMC.R`

`network_est` takes an interaction matrix and a phylogenetic tree as input, depending on the model chosen, and returns a list of input objects and the MCMC samples of the specified parameters interest.

## Usage

`network_est(Z, slices = 10, tree = NULL, model.type = c('full', 'distance', 'affinity'), uncertainty = FALSE, ... )`

## Arguments
+ `Z`:  an H x J binary or count matrix of interactions between two sets of species, H and J. If `Z` is count it will be converted to binary.
+ `slices`: the number of slices (MCMC samples) to run. In the case of the full model or the phylogeny-only model a single sample/iteration is the average of H samples from H conditional distributions, one for each row of `Z`; as in the ICM model. For the affinity-only model sampling from the full joint is possible; see 'Details' for more information.
+ `model.type`: either distance, affinity, or full to combine both distance and affinity. Default is the full model.
+ `uncertainty`: whether to sample and incorporate the uncertainty parameter or not. See 'Details' for more information.
+ `...`: optional arguments to the lower layer MCMC sampling algorithm:
  + `y`: a single number of a vector of size H specifying the initial value for the row affinities, default is 1. 
  + `a_y`: alpha hyperparameter of the Gamma priors of the row affinities, default is 1. Beta parameter is set to 1. 
  + `y_sd`: initial standard deviation for the proposal distribution (Gaussian) for the row parameters, default is 0.2. 
  + `w`: a single number of a vector of size J specifying the initial value for the column affinities, default is 1. 
  + `a_w`: alpha hyperparameter of the Gamma priors of the column affinities, default is 1. Beta parameter is set to 1.
  + `w_sd`: initial standard deviation for the proposal distribution (Gaussian) for the column parameters, default is 0.2.
  + `eta`: initial value for the tree transformation parameter under the `EB` model, default is 0.
  + `eta_sd`: initial standard deviation for the proposal distribution (Gaussian) for the column parameters, default is 0.005.
  + `burn.in`:  proportion of the initial samples to discard (0 to 1), default is 0.5.
  + `batch.size`: update frequency of the Adaptive MH for `_sd` parameters, default every 50 iterations. This applies to the affinity-only model.
  + `beta`:  tuning parameter for updating the Adaptive MH (>=0, 0 for no tuning, larger means more aggressive tuning).


## Details

The function `network_est()` returns a list that includes a list of posterior samples for each estimated parameter with burn-in removed. In addition, the input data; `Z` and `tree` for `full` and `distance` models, and `Z` only for `affinity` model are included in the output. Note that `network_est()` does initial cleaning of `Z` and `tree` to conform to the needed input. Such as:
    + converting `Z` to binary;
    + removing row-species in `Z` that are not in `tree`;
    + removing empty columns from `Z`;
    + removing tips in `tree` that do not correspond to row-species in `Z`;
    + ordering rows of Z to conform with transformation of `tree` to a distance matrix by `cophenetic`.

For a more details refer to `network_clean()` function. 


For convenience, row parameters (gammas), and column parameters (rhos), are indicated by the variables `y` and `w`, respectively, in the output. If using the `distance` model or `full`, this will also include posterior samples for the phylogenetic transformation parameter `eta`. The uncertainty parameter is denoted by `g`, and returned only when `uncertainty=TRUE`, otherwise `NULL`.

Specifying initial values for affinity parameters and related options in `...` is only relevant for the `full` or `affinity` model, otherwise they are ignored. Initial values for the `eta` parameter only apply in the `full` or `distance` model.
    
`network_est()` runs two types of upper layer MCMC sampling functions. It runs an MCMC using the full joint distribution when the `affinity` model is used, otherwise it runs an iterated conditional modes (ICM) MCMC sampler when the `full` or `distance` models are used. A mode of the ICM is the conditional joint distribution of a row on all other rows of `Z`. Besag(74, Sec 6.1) for details on the ICM method.

## References

Besag, J. (1974). _Spatial interaction and the statistical analysis of lattice systems._ Journal of the Royal Statistical Society. Series B (Methodological), 192–236.

Elmasri, M., Farrell, M., and Stephens, D. A. (2017). _A hierarchical Bayesian model for predicting ecological interactions using evolutionary relationships._ [arXiv](https://arxiv.org/abs/1707.08354).

Fritz, S. A., O. R. P. Bininda-Emonds, and A. Purvis (2009). Geographical variation in
predictors of mammalian extinction risk: big is bad, but only in the tropics. Ecology
letters 12 (6), 538–549.

Stephens, P. R., P. Pappalardo, S. Huang, J. E. Byers, M. J. Farrell, A. Gehman, R. R.
Ghai, S. E. Haas, B. Han, A. W. Park, J. P. Schmidt, S. Altizer, V. O. Ezenwa, and C. L.
Nunn (2017). _Global Mammal Parasite Database version 2.0._ Ecology 98 (February),
2017.

## Examples
A direct example from Elmasri et al. (2017) using the Global Mammal Parasite Database version 2.0 (GMPD).

```R
## Loading required packages
## Loading required packages
library(ape)
library(geiger)
library(fulltext)

## loading mammal supertree included in Fritz et al. (2009)
source('example-GMPD/download_tree.R')       # see variable 'tree'
## trimming 80% of tree tips for speed
pruned.tree <- drop.tip(tree,
                        sample(tree$tip.label)[1:(0.8*length(tree$tip.label))]) # trimming 80% of tree tips for speed

## loading GMPD
source('example-GMPD/load_GMPD.R')           # see matrix 'com'

## sourcing MCMC script
source('network_MCMC.R')

## running the model of interest
> obj = network_est(Z = com, slices=1000, tree=pruned.tree, model.type='full') # full model
Warning: Z is converted to binary!
Warning: not all species in tree exist in Z; missing are removed from tree!
Warning: not all row-species in Z exist tree; missing are removed from Z!
[1] "Ordering the rows of Z to match tree, and left ordering the columns.."
[1] "Running full model..."
[1] "Run for 1000 slices, and 500 burn-ins"
[1] "Slices 100"
[1] "Slices 200"
[1] "Slices 300"
[1] "Slices 400"
[1] "Slices 500"
[1] "Slices 600"
[1] "Slices 700"
[1] "Slices 800"
[1] "Slices 900"
[1] "Slices 1000"
[1] "Done!"
> names(obj)
[1] "param" "tree"  "Z"    
> names(obj$param)
[1] "w"       "y"       "eta"     "g"       "burn.in" "sd"  

## load useful network analysis functions
source('network_analysis.R')

## plot input matrix and tree
plot(obj$tree, show.tip.label=FALSE)
plot_Z(obj$Z)
```

Phylogenetic tree          |  Species interaction matrix Z
:-------------------------:|:-------------------------:
![](https://github.com/melmasri/HP-prediction/blob/master/example-GMPD/tree_example.png)  |  ![](https://github.com/melmasri/HP-prediction/blob/master/example-GMPD/Z_example.png)


```R
## Creating the H x J probability matrix
## Extracting mean posteriors
Y = if(is.matrix(obj$param$y)) rowMeans(obj$param$y) else  mean(obj$param$y)
W = if(is.matrix(obj$param$w)) rowMeans(obj$param$w) else  mean(obj$param$w)
Eta = if(is.null(obj$param$eta)) 0 else mean(obj$param$eta)

## Creating posterior distance
distance = 1/cophenetic(rescale(obj$tree, 'EB', Eta))
diag(distance)<-0
distance = distance %*% obj$Z
distance[distance==0]<-1

## Probability matrix
P = 1-  exp(-outer(Y, W)*distance)

## ROC curves, AUC and posterior Z
roc = rocCurves(obj$Z, obj$Z, P = P, all = TRUE) # ROC
plot_Z(P>roc$threshold + 0)                      # posterior Z

```

ROC curve and AUC     |  Posterior of Z
:-------------------------:|:-------------------------:
![](https://github.com/melmasri/HP-prediction/blob/master/example-GMPD/roc_example.png)  |  ![](https://github.com/melmasri/HP-prediction/blob/master/example-GMPD/Zpost_example.png)


To run the same analysis on the whole GMPD, use `tree` as is without pruning. For more information see `main_network.R`.


## Cross validation

For a cross validation example refer to scrip `mainCV_network.R`.
