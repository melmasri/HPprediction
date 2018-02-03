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
source('example-GMPD/')
com = readRDS('../com_oct25_2017.rds')
## sourcing MCMC script
library(microbenchmark)
source('network_MCMC.R')
## source('temp_fun.R')
## sourceCpp('network_MCMC.cpp')

cleaned = network_clean(com, tree, 'full', uncertainty=FALSE)
com = cleaned$Z                         # cleaned binary interaction matrix
tree = cleaned$tree                     # cleaned tree

param = network_est(com, slices=1000, tree=tree, 'aff', sparse=TRUE, backup=TRUE)

microbenchmark(network_est(com, slices=1, tree=tree, 'aff', sparse=TRUE), times=1)

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
[1] ICM_est(unname(com), tree, slices = 1, burn.in = 0) 444.5624 444.5624 444.5624
[2] ICM_est(unname(com), tree, slices = 1, burn.in = 0) 376.094 376.094 376.094
[3] ICM_est_C(unname(com), tree, slices = 1, burn.in  ) 63.05724
[4] ICM_est_C(unname(com), tree, slices = 1, burn.in )  163.1953
[5] ICM_est_C(unname(com), tree, slices = 1, burn.in)   344.6991
[6] ICM_est_C(unname(com), tree,  = 0, sparse = TRUE) 166.504

[1] old
[2] phangorn:::coph + no which
[3] 1 & 2 + sparse + only affinity
[4] 1 & 2 + sparse + affinity + dist (without assignment)
[5] 1 & 2 + sparse + affinity + dust with assignment
[6] 1 & 2 + sparse + no row extract + as matrix with assignment

Z = unname(com)
ind <- which(Z==1, arr.ind=TRUE)
sparseZ = sparseMatrix(ind[,1], j=ind[,2], x= rep(1, nrow(ind)),  dims=dim(Z))
dist = cophenetic(tree)
i=1
pdist.new = dist[i,]%*%Z
pdist = dist%*%Z

microbenchmark(pdist[1, ]<-pdist.new, dpdist[1,]<-dpdist.new,
               dpdist@x[ind + (i-1)]<-dpdist.new,
               dpdist@x[ind]<-dpdist.new@x)

microbenchmark(pdist[1, 1:10]<-1, dpdist[1,1:10]<-1, dpdist@x[ind + (i-1)]<-1,
               dpdist@x[ind]<-1, unit='ms', times=100)

microbenchmark(as(dpdist, 'matrix'), unit=' s')
i=1
dpdist.new = dist[i,]%*%sparseZ
dpdist = dist%*%sparseZ



Z = unname(com)
U0 = matrix(rexp(prod(dim(Z))), nrow(Z), ncol(Z))
Zs = Matrix(Z, sparse= TRUE)

microbenchmark(rg(Z, U0), rg(Zs, U0), times=10, unit='s')

rg<-function(Z,l){
    ## Sampling the uncertainty parameter
    ## Method one: P(z=0|g) = 1\delta_{s=0} + g\delta_{s>0}; dirac delta
    ZZ = 1*(l<1)                  # 1-inflated exponential
    ## Method 1
    M = sum(ZZ*Z)              # l<1(OR S>0) and Z=1,  N++
    N = sum((1-Z)*ZZ)          # l<1(OR S>0) and Z=0,  N-+
    rbeta(1 , N + 1, M + 1)         
}

n = 100
D = matrix(rExp(n^2), n,n)
A = diff(D)

microbenchmark(colMeans(diff(t(D))),rowMeans(D[,1:(n-1)] - D[,2:n] ))
               
library(Matrix)
library(microbenchmark)

Z = unname(com)
U0 = matrix(rexp(prod(dim(Z))), nrow(Z), ncol(Z))
Zs = Matrix(Z, sparse= TRUE)
Z0 = which(Z==0)
Z01 =   Z==0
unif = matrix(runif(length(U0)), nrow(U0), ncol(U0))

identical(rExp2(U0, 0.1, Z, Z01), rExp3(U0, 0.1, Z, Z0))

microbenchmark(rExp2(U0, 0.1, Z, Z01), rExp3(U0, 0.1, Z, Z0), times=50)

rExp2<-function(l, g, Z, Z0){
    ## Sampling from a zero-inflated Gumbel or a one-inflated Exponential
    ## Method one: P(z=0|g) = 1\delta_{s=0} + g\delta_{s>0}; dirac delta
    p = 1- exp(-l)
    ## unif = matrix(runif(length(l)), dim(p))
    U = 1 + 0*p
    U[!Z0] = -log(1- unif[!Z0]*p[!Z0])/(l[!Z0] + tol.err)
    aa = Z0 & (unif < g*p/(g*p + 1-p))
    U[aa] =  -log(1 - unif[aa]*(g*p[aa] + 1-p[aa])/g)/l[aa]
    U
}

