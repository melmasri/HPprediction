#########################################

## General variables
## please specify the following parameters
SAVE_PARAM = TRUE                    # should workspace be saved
SAVE_FILE = 'param.RData'            # name of output R workspace file
MODEL = 'full'                       # full, distance or affinity
SLICE = 100                          # no of iterations
subDir = ''                          # directory to print the results 
COUNT = TRUE                         # TRUE = count data, FALSE = year of first pub.

## loading mammal supertree included in Fritz et al. (2009)
source('example-GMPD/download_tree.R')       # see variable 'tree'

## loading GMPD
com = readRDS('../com_oct25_2017.rds')
## sourcing MCMC script
library(microbenchmark)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(BH)
library(Matrix)
source('network_MCMC.R')
source('temp_fun.R')
sourceCpp('network_MCMC.cpp')

cleaned = network_clean(com, tree, 'full', uncertainty=FALSE)
com = cleaned$Z                         # cleaned binary interaction matrix
tree = cleaned$tree                     # cleaned tree

param  = ICM_est(unname(com),tree,slices=1, burn.in=0)
param  = ICM_est_C(unname(com),tree,slices=1,burn.in=0, sparse= TRUE)


microbenchmark(ICM_est(unname(com),tree,slices=1, burn.in=0), times=1)

microbenchmark(ICM_est_C(unname(com),tree,slices=3,sparse=TRUE), ICM_est_C(unname(com),tree,slices=1,sparse=FALSE),times = 3)


## load useful network analysis functions
library(microbenchmark)
library(Rcpp)
library(RcppArmadillo)
library(BH)                             # for boost C++ libraries
sourceCpp('~/Github/HP-ICM/network_MCMC.cpp')


n = 1000
nrow = 300
w = rexp(n)
Z = matrix(rbinom(n*nrow, 1, prob=0.1), nrow, n)
Udp = matrix(rexp(nrow^2), nrow)%*%Z
hyper = c(0.1, 1)
i= 1
z = Z[i,]
Ud = Udp[i, ]
sig = rgamma(n, 0.1)

microbenchmark(raffinity.MH(w, z, Ud, sig, hyper),                         raffinity_MH_Rcpp(w, z, Ud, sig, hyper))


Unit: seconds
                                     expr      min       lq     mean   median
ICM_est_C(unname(com), tree, slices = 3) 170.4763 171.2151 171.7833 171.4158
ICM_est_C(unname(com), tree, slices = 3) 165.6226 167.2747 167.7588 167.6448
       uq      max neval
 171.5033 176.1518    10
 168.1898 170.3523    10





Code submitted to JASA with new dataset
Unit: seconds
                                                expr      min       lq     mean
ICM_est(unname(com), tree, slices = 1, burn.in = 0) 444.5624 444.5624 444.5624
above + phangorn 412.436 412.436 412.436

