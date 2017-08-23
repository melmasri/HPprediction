##Example for HP-Prediction (Elmasri et al. 2017)##

This example uses a small subset of the GMPD 2.0 ( DOI: 10.1002/ecy.1799) to show how to run the models described in Elmasri et al. 2017 (https://arxiv.org/abs/1707.08354).


### **Accessing data**

We use functions from the fulltext package (https://github.com/ropensci/fulltext) to directly download host-parasite interaction data and the host phylogeny. 

To download and format the GMPD 2.0 and subset to a manageable size for this example we run the following:

```r
source("load_GMPD.R")
com <- com[1:25,]
com <- com[,colSums(com)>0]
```
This results in a presence-absence community matrix with 
hosts as rows, and parasites as columns. 

For our paper we use the mammal supertree included in Fritz et al. 2009 (DOI: 10.1111/j.1461-0248.2009.01307.x).

To download a copy of the tree to your working directory we run the following script:

```r
source("download_tree.R")
```

### **Running the Models**

First we source the MCMC algorithms

```r
source('../networkMCMC.R')
```

All model variants can be run using the wrapper function network_est()

### Basic arguments ###

**Z** - The community matrix
**tree** - The phylo object
**model.type** - Either "distance", "affinity", or "both"
**slices** - This is the number of iterations to run the model
**uncertainty** - Whether to extend model to estimate uncertainty in unobserved interactions (TRUE/FALSE), default is FALSE 

You can pass the following optional arguments:

### Priors ###

**a_w** - the alpha hyper parameter of the Gamma distribution for parasite parameters

**a_y** - the alpha hyper parameter of the Gamma distribution for host parameters

**w_sd** - the initial variance for the proposal distribution (Gaussian) for the host parameters

**y_sd** - the initial variance for the proposal distribution (Gaussian) for the host parameters


### Burn in set-up ###

**burn.in** - Proportion of the initial samples to discard (0 to 1), default is 0.5

**batch.size** - The update frequency of the Adaptive MH for _sd parameters, default every 50 iterations. This applies to affinity-only model.

**beta** - Tuning parameter for updating the Adaptive MH ( >=0, 0 is no tuning, larger means more aggressive tuning).


### Additional arguments for the full or distance-only model ###

In addition to the ones above you have

**eta** - The initial starting value of eta, default is 0

**eta_sd** - Initial variance for the proposal distribution (Gaussian) for eta, default is 0.005


### Full model ###
To run the full model (layering of affinity and phylogeny models) without the uncertainity extension, we use the following code: 

```{r
obj <- network_est(Z = com, slices=100, tree=tree, model.type='both')
```

This function returns a list that includes the posterior for each estimated parameter, with the burn-in removed.

For convenience, host parameters (gammas), and parasite parameters (rhos), are indicated by the variables w and y in the output. If using the distance-onnly or full method, this will also include the posterior for the phylogenetic transformation parameter eta. 

To convert these estimates to interaction probabilities we first take the posterior means for the host, parasite, and eta parameters, and plug them into equation 4 in Elmasri et al. 2017.

```{r
Y <- if(is.matrix(obj$param$y)) rowMeans(obj$param$y) else  mean(obj$param$y)

W <- if(is.matrix(obj$param$w)) rowMeans(obj$param$w) else  mean(obj$param$w)

Eta <- if(is.null(obj$param$eta)) 0 else mean(obj$param$eta)

distance = 1/cophenetic(rescale(obj$tree, 'EB', Eta))
diag(distance)<-0
distance = distance %*% obj$Z
distance[distance==0]<-1

P = 1 - exp(-outer(Y, W)*distance)

```

P is the interaction probability matrix for the full model. 


For the affinity only model the P matrix is calculated as:

```r
P = 1 - exp(-outer(Y, W))
```

For the distance only model the P matrix is calculated as:
```r
P = 1 - exp(-distance)
```
