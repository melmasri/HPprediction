rm(list=ls())
source('library.R')
source('gen.R')

## Global Variable
DATAFILENAME = 'comGMPD.RData'
print(DATAFILENAME)
source('library.R')
source('gen.R')
load(DATAFILENAME)

#######################
## subsetting
    ##GMP
aux = rownames(com) %in% pan$bionomial[pan$Order=="Carnivora"]
phy_dist = phy_dist[aux,]
phy_dist = phy_dist[,aux]
diag(phy_dist)<-0
rownames(phy_dist)<-colnames(phy_dist)<-NULL
range(phy_dist[upper.tri(phy_dist)])


D = generate_interactions(70, 500, eta = 1.5, aj=2, ai=0.7, D = phy_dist)
plot_Z(D$Z)
dim(D$Z)
hyper = list(parasiteHyper =c(2, 1), hostHyper = c(0.7,1), etaSamplingSD = c(0.01)) 
param = ICM_est(D$Z,slice=10000,hyper = hyper, dist = D$D^D$eta,AdaptiveMC=TRUE,ICM.HORIZ= TRUE, updateHyper = FALSE)
ana.plot(param)
aux = getMean(param)
r = which.max(rowMeans(param$y));D$y[r]
c = which.max(rowMeans(param$w));D$w[c]
aux$y[r];aux$w[c];aux$y[r]*aux$w[c]
rowMeans(param$hh)
plot(param$w[c,]*param$y[r,])
abline(h=D$w[c]*D$y[r])
plot(param$hh[1,])
aux$eta=D$eta
P = P = 1-  exp(-outer(aux$y, aux$w))
P = 1-  exp(-outer(aux$y, aux$w)*((D$D^aux$eta)%*%D$Z))

P = 1-  exp(-outer(D$y, D$w)*((D$D^D$eta)%*%D$Z))

roc.all = rocCurves(Z=D$Z, Z_cross= D$Z, P=P, plot=FALSE, bins=400, all=TRUE)
tb.all  = ana.table(D$Z,D$Z, roc=roc.all, plot=FALSE)
tb.all




## Converting to Matlab
DATAFILENAME = 'comGMPD.RData'
load(DATAFILENAME)
library(R.matlab)


Z=1*(com>0)
## Testing S
source('gen.R')
source('library.R')

hyper = list(parasiteHyper =c(0.32, 1), hostHyper = c(0.94,1), etaSamplingSD = c(0.01))

slice = 1000
ICM.HORIZ = TRUE
param = ICM_est(Z =1*(com>0),slice=slice ,dist= phy_dist,
    hyper = hyper, AdaptiveMC=TRUE,ICM.HORIZ= ICM.HORIZ, eta=1)
aux = getMean(param)

P = 1-  exp(-outer(aux$y, aux$w)*exp(-aux$eta* (phy_dist)%*%(-log(param$U))))
P = 1-  exp(-outer(aux$y, aux$w)*((phy_dist^aux$eta)%*%Z))
P = 1-  exp(-outer(aux$y, aux$w)*((phy_dist^aux$eta)%*%(-log(param$U))))
range(P)
roc.all = rocCurves(Z=Z, Z_cross= Z, P=P, plot=TRUE, bins=400, all=TRUE)
tb.all  = ana.table(Z,Z, roc=roc.all, plot=TRUE)
tb.all

 exp with no Eta
    auc    thresh tot.inter hold.out pred pred.all
1 72.76 0.1629073      3966        0  NaN 0.704236

Using Exp with Eta
auc    thresh tot.inter hold.out pred  pred.all
1 71.47 0.1904762      3966        0  NaN 0.7748361
> ana.plot(param)
> 


Using Z with no ETA
 auc     thresh tot.inter hold.out pred  pred.all
1 89.6 0.03258145      3966        0  NaN 0.7965204

Using Z with Eta, Eta est = 1.24
  auc     thresh tot.inter hold.out pred  pred.all
1 89.8 0.02255639      3966        0  NaN 0.8055976

Using S with no Eta
  auc     thresh tot.inter hold.out pred  pred.all
1 86.92 0.03508772      3966        0  NaN 0.7428139

Using S with Eta 
   auc     thresh tot.inter hold.out pred  pred.all
1 86.86 0.01002506      3966        0  NaN 0.7854261
> 




#### Hamiltonian MC
n=1000
x = cbind(rnorm(n, 0,1) ,rnorm(n, 50,1), rnorm(n, 3,1))
U <-function(q)    sum((t(x)-q)^2)/2 + sum(q^2)/2
grad_U <-function(q) -rowSums(t(x)-q) + q

vec = matrix(1, ncol(x), 10000)
for(i in 2:ncol(vec))
    vec[,i] = HMC(U, grad_U, epsilon = 0.001, L=10, vec[,i-1])
vec[,ncol(vec)]
plot(t(vec[,-c(1:500)]))





x = rnorm(1000,5,1)
U <-function(q)    sum((x-q)^2)/2 + sum(q^2)/2
grad_U <-function(q) -sum(x-q) + q

vec=rep(0,10000)
for(i in 2:length(vec))
    vec[i] = HMC(U, grad_U, epsilon = 0.01, L=10, vec[i-1])
vec[10000]
plot(vec)

log.posterior(c(y0in[i, 1], petain[i], w0in[,i]),d=dist[i,], Z=Z, z=Z[i,],
              nh=mr[i],U=U0[i,], a_y=hh[1,s], a_w=hh[3,s])
grad.log.posterior(c(y0in[i, 1], petain[i], w0in[,i]),d=dist[i,], Z=Z, z=Z[i,],
                   nh=mr[i],U=U0[i,], a_y=hh[1,s], a_w=hh[3,s])

HMC(U=log.posterior, grad_U = grad.log.posterior, epsilon=0.001,L=15,c(y0in[i, 1], petain[i], w0in[,i]),d=dist[i,], Z=Z, z=Z[i,],nh=mr[i],u=U0[i,], a_y=hh[1,s], a_w=hh[3,s] )

log.posterior<-function(q, ...)
{
    el<-list(...)
    eta=q[2]
    y=q[1]
    w = q[-c(1,2)]
    ## Condtional row posteriors
    pdist=(el$d^eta)%*%el$Z + tol.err
    posterior =  el$nh*log(y) + sum(el$z*(log(w) + log(pdist))) - sum(el$u*y*w*pdist)
    prior = (el$a_y-1)*log(y) - y +
        (el$a_w-1)*sum(log(w)) - sum(w) -log(eta)
    -(posterior + prior)
}

grad.log.posterior<-function(q, ...)
{
    el<-list(...)
    eta=q[2]
    y=q[1]
    w = q[-c(1,2)]
    pdist=(el$d^eta)%*%el$Z + tol.err
    pU = pdist*el$u
    dy = -el$nh/y + sum(w*pU) - (el$a_y-1)/y + 1
    dw = -el$z/w + y*pU - (el$a_w-1)/w +1
    deta = log(eta)*(-el$nh+ sum(y*w*pU)) +1/eta
    c(dy, deta, dw)
}


log.posterior(c(y0in[i, 1], petain[i], w0in[,i]),d=dist[i,], Z=Z, z=Z[i,],
              nh=mr[i],U=U, a_y=  hh[1,s], a_w = hh[3,s])

new.eta = rEta.Horiz(petain[i],dist,pdist[i,],i,Z,y0in[i,1],
                            w0in[,i+1],U0[i,], eta_sd=eta_sd)

test<-function(a, b,...){
    el <-match.call()
    print(el)
    a*b*el$c
}
x =1;z=3;c=4
test(1,2,x,z,c=5)

raffinity.MH<-function(old, z, Ud, sig=0.1, hyper){
    a = hyper[1]
    b = hyper[2]
    epsilon = rnorm(length(old), 0, sd = sig)
    new = old +  -2*(old + epsilon <=0)*epsilon + epsilon    

    likeli = (z+a-1)*log(new/old) -  (new - old)*(b +Ud)
    
    ratio = pmin(1, exp(likeli))
    u = 1*(runif(length(old))<=ratio)
    u*new + (1-u)*old
}


Z = com_paCross
slice=slice
dist= dist
eta=1,
length(param$eta)
aux$eta = param$eta[480]
aux$w = param$w[,480]
aux$y = param$y[,480]

DATAFILENAME = 'comGMPD.RData'
load(DATAFILENAME)
source('library.R')
cnames = colnames(com)
rnames = rownames(com)
com = unname(com)
phy_dist = unname(phy_dist)

pairs = cross.validate.fold(1*(com>0), n=5)
com_pa=  1*(com>0)
Z = com_pa
x=1
ICM.HORIZ= TRUE
com_paCross = Z
com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
hyper = list(parasiteHyper =c(0.32, 1), hostHyper = c(0.94,1), etaSamplingSD = c(0.01)) 
which(phy_dist%*%com_paCross==0,arr.ind=TRUE)

source('HM.R')
source('gen.HM.R')
slice = 2000
dist = phy_dist
diag(dist)=0

param = ICM_est(com_paCross,slice=slice ,dist= dist,eta=0.1,
    hyper = hyper, AdaptiveMC=FALSE,ICM.HORIZ= ICM.HORIZ)

aux = getMean(param);aux
ana.plot(param)

testfunc<-function(q, ...)
    log.posterior.full(q, ...)
P =1-exp(-outer(aux$y, aux$w)*(1 + aux$eta*dd/Zmc))
dd = exp(-aux$eta*(dist%*%com_paCross))
dd = (dist%*%com_paCross)/aux$eta
dd = (dist%*%com_paCross)
P = 1-  exp(-outer(aux$y, aux$w)*dd)
range(P)
roc = rocCurves(Z=com_pa, Z_cross= com_paCross, P=P, plot=TRUE, bins=400, all=FALSE)
tb  = ana.table(com_pa, com_paCross, roc=roc, plot=FALSE)
tb
P = 1-exp(-phy_dist%*%com_paCross)
P = 1-exp(-d1%*%com_paCross)
P = 1-exp(-(com_paCross%*%t(com_paCross))%*%com_paCross)
param$y=param$y[-c(1:400),]
param$w=param$w[-c(1:400),]
param$eta = param$eta[-c(1:400)]

## Standard w=1;y=1;eta=1;pdist=dist%*%Z/eta
    auc    thresh tot.inter hold.out      pred pred.all
1 87.21 0.4310777      3966      694 0.7752161 0.677761

### After fixing HMC - affinity only, without phi_hj not dist at all
tb
    auc     thresh tot.inter hold.out      pred  pred.all
1 77.61 0.01503759      3966      692 0.6748555 0.6913767
## affinity, with dist-without sim
> 
 tb
   auc     thresh tot.inter hold.out      pred  pred.all
1 89.5 0.01503759      3966      692 0.8439306 0.8050933



### Stan testign
library(rstan)
model <- stan_model(file="stan/HP_subIter.stan")

stan.data <- list(H=n_y, P=n_w, z=Z[i,], n=mr[i], u=U0[i,], alpha_h=hh[1,s], alpha_p=hh[2,s])

iter=5
fit <- sampling(model, data=stan.data,iter=iter, chains=1, cores=1, thin=iter,check_data= FALSE)
extract(fit, include = TRUE)$q

system.time(
    sapply(1:iter, function(r) HMC(U=log.posterior.affinity, grad_U = grad.log.posterior.affinity,
                    epsilon=0.01,L=10,
                    c(log(y0[i, s]),log(w0[,s])),
                    d=dist[i,], Z=Z, z=Z[i,],nh=mr[i],u=U0[i,],
                    a_y=hh[1,s], a_w=hh[2,s], dd=dd[i,], pdist=1)))





a = runif(5)
b= runif(8)

a%*%t(b)



### testing no empty columns in grouping
com_paCross = com_pa
x=1
com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0

sapply(1:tot.gr, function(i){
    A = 1*(com_pa>0)
    A[pairs[which(pairs[,'gr']==i),c('row', 'col')]]<-0
    which(phy_dist%*%A==0, arr.ind=TRUE)
})

a = unique(aux[,'col'])
colSums(A[,a])
which(colSums(A)==1)

which(colSums(A)==0)


cbind(phy_dist[188, ],A[, 208])



library(ape)


dist = phy_dist
dd = 1/dist
diag(dd)<-0

                                        #diag(dd)<-max(dd)
Z = com_pa
dd= (dd/max(dd) - 1)%*%Z


load('comEID-
    PS.RData')
Z= 1*(com>0)
pairs = cross.validate.fold(1*(com>0), n=5)
com_paCross = 1*(com>0)
com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
Z=com_paCross

grid=seq(0.1,1,0.05)
aux =sapply(grid, function(eta){
    print(eta)
    ZtZ = Z%*%t(Z)
    phy_dist=1/phy_dist
    
    diag(ZtZ)<-0
    pdist = ZtZ%*% Z
    P = 1-exp(-pdist)
    roc = rocCurves(Z=1*(com>0), Z_cross= com_paCross,
        P=P, plot=FALSE, bins=400, all=FALSE)
    tb  = ana.table(1*(com>0), Z, roc=roc, plot=FALSE)
    cbind(eta=eta, tb=tb)
})

png('lambda_trans_GMPD-geiger.png')
plot(unlist(aux['eta',]), unlist(aux['tb.auc',]), ylab = 'AUC', xlab='parameter')
dev.off()

aux[,which.max(unlist(aux['tb.auc',]))]


load('comGMPD.RData')
tree <- read.tree('../Data/mammals.tre')
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(com)])

# Testing all names in hosts in com exist in tree
if(! all(sapply(rownames(com), function(r) r %in% tree$tip.label))) {
		print('Warning! Not all hosts in com exist in tree. Hosts not found in tree will be removed.')
		com <- com[rownames(com)%in%tree$tip.label,]
			}

# Scaling the tree.
phy_dist<- cophenetic(rescale(tree, "a", 2))
phy_dist = dist_ordering(phy_dist, com)

phy_dist<- cophenetic(rescale(tree, "kappa", 1))

phy_dist<- cophenetic(rescale(tree, "lambda", 1))



library(geiger)
source('library.R')
load('comGMPD.single.RData')
tree <- read.tree('../Data/mammals.tre')
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(com)])

grid=seq(-.05,0.05,0.001)
Z= 1*(com>0)
aux =sapply(grid, function(eta){
    print(eta)
    phy_dist<- cophenetic(rescale(tree, "EB", eta))
    phy_dist = dist_ordering(phy_dist, com)
    dd =1/phy_dist
    diag(dd)<-0
    pdist = dd %*% Z
    P = 1-exp(-pdist)
    roc = rocCurves(Z=Z, Z_cross= Z, P=P, plot=FALSE, bins=400, all=TRUE)
    tb  = ana.table(Z, Z, roc=roc, plot=FALSE)
    cbind(eta=eta, tb=tb)
})

plot(unlist(aux['eta',]), unlist(aux['tb.auc',]), ylab = 'AUC', xlab='parameter')


aux[,which.max(unlist(aux['tb.auc',]))]

roc = rocCurves(Z=1*(com>0), Z_cross= Z, P=P, plot=TRUE, bins=400, all=FALSE)
tb  = ana.table(1*(com>0), Z, roc=roc, plot=FALSE)




## 1 + eta*tree
library(geiger)
source('library.R')
load('comGMPD.RData')
tree <- read.tree('../Data/mammals.tre')
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(com)])
phy_dist<- cophenetic(rescale(tree, "lambda", 1))
phy_dist = dist_ordering(phy_dist, com)
#phy_dist = phy_dist/max(phy_dist)
phy_dist= 1/phy_dist
diag(phy_dist)<-0
Z=1*(com>0)
grid=seq(0,10,0.1)

aux =sapply(grid, function(eta){
    print(eta)
    dd = (phy_dist^eta)%*%Z
    pdist = dd
    P = 1-exp(-pdist)
    roc = rocCurves(Z=Z, Z_cross= Z, P=P, plot=FALSE, bins=400, all=TRUE)
    tb  = ana.table(Z, Z, roc=roc, plot=FALSE)
    cbind(eta=eta, tb=tb)
})
aux



##
Z= 1*(com>0)
Z = unname(Z)
phy_dist = unname(phy_dist)
colSums(Z)
col = 743
Z[241, col]<-0

eta = 2
dist = phy_dist
pdist = (dist^eta)%*%Z
pdist0 = pdist==0

which(pdist0, arr.ind=TRUE)

pdist.zero= pdist0
alpha.old=0.5
alpha_sd= 0.1
y = rep(1,nrow(Z))
w = rep(1, ncol(Z))
U = rExp(pdist*outer(y,w))
U[Z==0]<-1
    
rAlpha<-function(alpha.old,n0,y, w, U, sd =0.01){
    ## MH for alpha, the minimum jump when there is no information
    alpha.prop = alpha.old*exp(rnorm(1, 0, sd = alpha_sd))
    alpha.prop =min(alpha.prop,5)
    
    likeli = n0*(log(alpha.prop) - log(alpha.old)) - 
        sum(outer(y,w)[pdist.zero]*U[pdist.zero])*(alpha.prop - alpha.old)
    
    u = (runif(1)<= min(1, exp(likeli)))
    if(u) alpha.prop else alpha.old
}


grid = seq(0,10, 0.1)
aux= sapply(grid, function(alpha.prop){
    likeli = n*(log(alpha.prop) - log(alpha.old)) - 
        sum(outer(y,w)[pdist.zero]*U[pdist.zero])*(alpha.prop - alpha.old)
   })
plot(grid, aux)
cbind(aux, grid)


library(geiger)

tree <- sim.bdtree(n=5)

plot(tree, edge.width = 2)

cbind(tree$edge, tree$edge.length)
tree2= rescale(tree, "EB", 1)
cbind(tree2$edge, tree2$edge.length)
cophenetic(tree2)

tree3= rescale(tree, "EB", -1)
cbind(tree3$edge, tree3$edge.length)
cophenetic(tree3)

> cbind(tree2$edge, tree2$edge.length)
     [,1] [,2]      [,3]
[1,]    6    7 1.5424207
[2,]    7    1 1.3095209
[3,]    7    8 0.3190863
[4,]    8    2 0.8615852
[5,]    8    9 3.2949152
[6,]    9    3 1.9040293
[7,]    9    4 3.8992498
[8,]    6    5 1.4450434
> 

cophenetic(tree)

tree2$edge[14:15,]
tree2$edge.length[14:15]

tree$edge[14:15,]
tree$edge.length[14:15]*exp(r)
names(tree)
tree
log(tree2$edge.length[14]/tree$edge.length[14])/r
which.min(tree$edge.length)

exp(r*tree$edge.length[2])
tree2$edge.length[2]
tree$edge.length[2]

tree$edge.length[1]
log(tree2$edge.length[1]/tree$edge.length[1])
log(tree2$edge.length[2]/tree$edge.length[2])
log(tree2$edge.length[3]/tree$edge.length[3])
log(tree2$edge.length[4]/tree$edge.length[4])



cophenetic(phy)

a= 0.0001

pdf('T_example_EB-1.pdf')
plot(rescale(phy, 'EB',-1))
dev.off()


cophenetic(
cophenetic(rescale(phy, 'EB', 0))

plot(rescale(phy, 'EB', 2))

sum(phy$edge.length[c(8,1,2 )])

exp(phy$edge.length[8]) + exp(sum(phy$edge.length[c(1,2)]))

## t1 - t4
(exp(a*sum(phy$edge.length[c(1,2)]))- exp(a*sum(phy$edge.length[c(1)])) +
exp(a*sum(phy$edge.length[c(1,3,4)])) - exp(a*sum(phy$edge.length[c(1)])))/a


cbind(phy$edge, phy$edge.length)
phy1 = rescale(phy, 'EB', 1)
cbind(phy1$edge, phy1$edge.length)

exp(sum(phy$edge.length[c(2)])) + 
exp(sum(phy$edge.length[c(3,4)])) 

phy = tree

phy1 = eb.phylo(phy)

eb.phylo=function(phy){


	ht=heights.phylo(phy)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), phy$edge[,2])
	ht$t1=Tmax-ht$end[phy$edge[mm,1]]
	ht$t2=ht$start-ht$end+ht$t1

	z=function(a, sigsq=1){
		if(a==0) return(phy)
		bl = (exp(a*ht$t2)-exp(a*ht$t1))/(a)
		phy$edge.length=bl[phy$edge[,2]]
        phy$edge.length=phy$edge.length*sigsq
		phy
	}
    
	attr(z,"argn")=c("a", "sigsq")
	return(z)
}

heights.phylo=function(x){
    ## from https://github.com/mwpennell/geiger-v2/blob/master/R/utilities-phylo.R#L1353-L1382
    
    phy=x
	phy <- reorder(phy, "postorder")
	n <- length(phy$tip.label)
	n.node <- phy$Nnode
	xx <- numeric(n + n.node)
	for (i in nrow(phy$edge):1) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
	root = ifelse(is.null(phy$root.edge), 0, phy$root.edge)
	labs = c(phy$tip.label, phy$node.label)
	depth = max(xx)
	tt = depth - xx
	idx = 1:length(tt)
	dd = phy$edge.length[idx]
	mm = match(idx, c(phy$edge[, 2], Ntip(phy) + 1))
	dd = c(phy$edge.length, root)[mm]
	res = cbind(start = tt+dd, end = tt)
	rownames(res) = idx
    data.frame(res)
}


 phy$edge
     [,1] [,2]
[1,]    9    3
[2,]    9    4
[3,]    8    2
[4,]    8    9
[5,]    7    1
[6,]    7    8
[7,]    6    7
[8,]    6    5
> nrow(phy$edge)



set.seed(100)

tree <- sim.bdtree(n=5,  seed=0)


cophenetic(rescale(rescale(tree, 'EB', 0), "depth", depth=1))
cophenetic(rescale(rescale(tree, 'EB', 0), "depth", depth=3))



plot(tree, edge.width = 2, main= 0)
add.scale.bar()

r = 5
tree2= rescale(tree, "lambda", 0.5)
cophenetic(tree2)
dev.new()
plot(tree2, edge.width = 2, main = r)

r = -5
tree3= rescale(tree, "EB", r)
cophenetic(tree3)
dev.new()
plot(tree3, edge.width = 2, main = r)



arrange.tree<-function(phy){
    ## taken from .b.phylo from
    ## https://github.com/mwpennell/geiger-v2/blob/master/R/utilities-phylo.R#L1353-L1382
    ht=heights.phylo(tree)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), phy$edge[,2])
	ht$t1=Tmax-ht$end[phy$edge[mm,1]]
	ht$t2=ht$start-ht$end+ht$t1
    ht
}

eb.phylo=function(phy, heights, a){
    ## modified form .eb.phylo
    ## https://github.com/mwpennell/geiger-v2/blob/master/R/utilities-phylo.R
	## #ht=heights.phylo(phy)
	## N=Ntip(phy)
	## Tmax=ht$start[N+1]
	## mm=match(1:nrow(ht), phy$edge[,2])
	## ht$t1=Tmax-ht$end[phy$edge[mm,1]]
	## ht$t2=ht$start-ht$end+ht$t1
    if(a==0) return(phy)
    bl = (exp(a*heights$t2)-exp(a*heights$t1))/(a)
    phy$edge.length=bl[phy$edge[,2]]
    phy
}


heights.phylo=function(x){
    phy=x
	phy <- reorder(phy, "postorder")
	n <- length(phy$tip.label)
	n.node <- phy$Nnode
	xx <- numeric(n + n.node)
	for (i in nrow(phy$edge):1) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
	root = ifelse(is.null(phy$root.edge), 0, phy$root.edge)
	labs = c(phy$tip.label, phy$node.label)
	depth = max(xx)
	tt = depth - xx
	idx = 1:length(tt)
	dd = phy$edge.length[idx]
	mm = match(1:length(tt), c(phy$edge[, 2], Ntip(phy) + 1))
	dd = c(phy$edge.length, root)[mm]
	ss = tt + dd
	res = cbind(ss, tt)
	rownames(res) = idx
	colnames(res) = c("start", "end")
	res = data.frame(res)
	res
}



library(geiger)
source('library.R')
load('comGMPD.RData')
tree <- read.tree('../Data/mammals.tre')
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(com)])

grid=seq(-.05,0.15,0.01)
Z= 1*(com>0)

aux =sapply(grid, function(eta){
    print(eta)
    phy_dist<- cophenetic(rescale(rescale(tree, "EB", eta), "depth", depth=1))
    phy_dist = dist_ordering(phy_dist, com)
    dd =phy_dist
    diag(dd)<-0
    pdist = dd %*% Z
    P = 1-exp(-pdist)
    roc = rocCurves(Z=Z, Z_cross= Z, P=P, plot=FALSE, bins=400, all=TRUE)
    tb  = ana.table(Z, Z, roc=roc, plot=FALSE)
    cbind(eta=eta, tb=tb)
})

plot(unlist(aux['eta',]), unlist(aux['tb.auc',]), ylab = 'AUC', xlab='parameter')

aux[,which.max(aux[2,])]



library(ape)
library(geiger)
load('~/Github/ICM/comGMPD.RData')
tree <- read.tree('~/Github/Data/mammals.tre')
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(com)])

# Testing all names in hosts in com exist in tree
if(! all(sapply(rownames(com), function(r) r %in% tree$tip.label))) {
		print('Warning! Not all hosts in com exist in tree. Hosts not found in tree will be removed.')
		com <- com[rownames(com)%in%tree$tip.label,]
}

source('library.R')
source('gen-pseudo.R')


dd = cophenetic(rescale(tree, 'EB', 0))
aux <- sapply(rownames(dd), function(r) which(r==rownames(com)))
com = com[aux,]
Z= unname(com)
Z = 1*(Z>0)

pairs = cross.validate.fold(1*(com>0), n=5,2)
Z[pairs[which(pairs[,'gr']==2),c('row', 'col')]]<-0

p = sqrt(-log(1- mean(rowMeans(1*(com>0)))))
a_w =a_y = p


param = gibbs_one(Z=Z,slice=5000,burn=0.5, a_w =p, a_y =p, eta = -.01 , tree = tree, distOnly= FALSE, eta_sd= 0.001)

aux = getMean(param);aux$eta
plot(rescale(tree, 'EB',aux$eta))
tree.ht = arrange.tree(tree)

yw = outer(aux$y, aux$w)
pdist= 1/cophenetic(rescale(tree, 'EB', aux$eta))
diag(pdist)<-0
pdist = pdist%*%Z
P = 1- exp(-yw*pdist)
roc = rocCurves(Z=1*(com>0), Z_cross= Z, P=P, plot=TRUE, bins=400, all=FALSE)
tb  = ana.table(1*(com>0), Z, roc=roc, plot=FALSE)
tb

plot_Z(1*(P>roc$thre))
aux$eta


> tb
    auc     thresh tot.inter hold.out pred  pred.all
1 87.81 0.02255639      3966        0  NaN 0.7665154
> plot_Z(1*(P>roc$thre))

## ay = 0.2 aw = 0.2
 tb
    auc     thresh tot.inter hold.out pred  pred.all
1 89.02 0.02255639      3966        0  NaN 0.8484619

## ay  0.4, aw = 0.2
 auc     thresh tot.inter hold.out pred  pred.all
1 88.9 0.02005013      3966        0  NaN 0.8497226

 tb  = ana.table(1*(com>0), Z, roc=roc, plot=FALSE)
> tb
   auc     thresh tot.inter hold.out pred  pred.all
1 88.7 0.02005013      3966        0  NaN 0.8386283

 ICM_est(Z=Z,slice=500,eta=0,tree=tree, burn=0.2,eta_sd = 0.01, a_w =0.1, a_y =0.1, beta=3)
 tb
    auc     thresh tot.inter hold.out pred pred.all
1 89.54 0.03258145      3966        0  NaN 0.818709
> 

param = ICM_est(Z=Z,slice=500,eta=0,tree=tree, burn=0.2,eta_sd = 0.01, a_w =0.01, a_y =0.01, beta=3)

 tb
   auc     thresh tot.inter hold.out pred  pred.all
1 89.7 0.03258145      3966        0  NaN 0.8408976

param = ICM_est(Z=Z,slice=500,eta=0,tree=tree, burn=0.2,eta_sd = 0.01, a_w =0.01, a_y =0.01, beta=5)
tb
    auc     thresh tot.inter hold.out      pred pred.all
1 88.38 0.02005013      3966      620 0.8596774 0.829299
> 

param = ICM_est(Z=Z,slice=500,eta=0,tree=tree, burn=0.2,eta_sd = 0.01, a_w =0.1, a_y =0.5, beta=5)

 tb  = ana.table(1*(com>0), Z, roc=roc, plot=FALSE)
> tb
    auc     thresh tot.inter hold.out      pred  pred.all
1 87.36 0.01253133      3966      620 0.8403226 0.8113969
> ana.plot(param)


param = ICM_est(Z=Z,slice=500,eta=0,tree=tree, burn=0.2,eta_sd = 0.01, a_w =0.1, a_y =0., beta=5)

 tb
    auc     thresh tot.inter hold.out     pred  pred.all
1 86.86 0.01002506      3966      620 0.866129 0.8330812

param = ICM_est(Z=Z,slice=500,eta=0,tree=tree, burn=0.2,eta_sd = 0.01, a_w =0.001, a_y =0.001, beta=5)

   auc     thresh tot.inter hold.out      pred  pred.all
1 88.55 0.02756892      3966      620 0.8419355 0.8146747
>

0.01
sqrt(mean(rowMeans(1*(com>0))))


-log(1- mean(rowMeans(1*(com>0))))
param = ICM_est(Z=Z,slice=500,eta=0,tree=tree, burn=0.2,eta_sd = 0.01, a_w =0.1, a_y =0.4, beta=5)
 auc     thresh tot.inter hold.out      pred  pred.all
1 87.92 0.02005013      3966      620 0.8112903 0.7866868

param = ICM_est(Z=Z,slice=500,eta=0,tree=tree, burn=0.2,eta_sd = 0.01, a_w =0.06, a_y =0.4, beta=5)
 auc     thresh tot.inter hold.out pred  pred.all
1 87.5 0.01253133      3966      620 0.85 0.8199697


param = ICM_est(Z=Z,slice=500,eta=0,tree=tree, burn=0.2,eta_sd = 0.01, a_w =0.4, a_y =0.06, beta=5)

> plot_Z(1*(P>roc$thre))

nrow(Z)/ncol(Z)

p = sqrt(-log(1- mean(rowMeans(1*(com>0)))))
p*nrow(Z)/ncol(Z)
p*ncol(Z)/nrow(Z)



 if(!uncertain){
                U0<- rExp(pdist*outer(y0[,i],w0[,i]))
                U0[Z0]<-1
            }else
                U0 <-rExp2(pdist*outer(y0[,i],w0[,i]), g0[i], Z, Z0)


if(uncertain){
    g0[i+1] = rg(Z, U0,l=pdist*outer(y0[,i],w0[,i])) 
}




OptionParser(option_list="-f comGMPD-year.single.RData -r uncertain --no.cores=1 --no.cycles=5000")

opt$runtype = 'uncertain'
opt$no.cycles = 5000
opt$no.cores=1
opt$help=FALSE
opt$file = 'comGMPD-year.single.RData'





load('~/Github/ICM/comGMPD-year.RData')
aux = rownames(com) %in% pan$bionomial[pan$Order=="Carnivora"]

com =com[aux,]
aux = colSums(1*(com>0))
com = com[,aux>1]
aux = rowSums(1*(com>0))
com = com[aux>1,]
com = lof(com)


plot_Z(com>0 + 0)


## Loading the tree
tree <- read.tree('../mammals.tre')
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(com)])

# Testing all names in hosts in com exist in tree
if(!all(sapply(rownames(com), function(r) r %in% tree$tip.label))) {
    print('Warning! Not all hosts in com exist in tree. Hosts not found in tree will be removed.')
    com <- com[rownames(com)%in%tree$tip.label,]
}


pdf('tree-GMPD.pdf')
plot(rescale(tree, 'EB', 0), font=1, cex=0.4, label.offset = 1, show.tip.label=FALSE)
dev.off()

pdf('tree_Carnivora-02.pdf')
plot(rescale(tree, 'EB', 0.02), font=1, cex=0.4, label.offset = 1, show.tip.label=FALSE)
dev.off()


pdf('tree_Carnivora-N02.pdf')
plot(rescale(tree, 'EB', -.02), font=1, cex=0.4, label.offset = 1, show.tip.label=FALSE)
dev.off()



load('~/Github/ICM/comGMPD.RData')


## Loading the tree
tree <- read.tree('../mammals.tre')
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(com)])

# Testing all names in hosts in com exist in tree
if(!all(sapply(rownames(com), function(r) r %in% tree$tip.label))) {
    print('Warning! Not all hosts in com exist in tree. Hosts not found in tree will be removed.')
    com <- com[rownames(com)%in%tree$tip.label,]
}



dd = cophenetic(rescale(tree, 'EB', 0))
host.order <- sapply(rownames(dd), function(r) which(r==rownames(com)))
com = com[host.order,]

pdf('Z_gmp.pdf')
plot_Z(lof(1*(com>0)), xlab = 'parasites', ylab='hosts')
dev.off()




grid1 = seq(3, 5, 0.01)

aa = sapply(grid1, function(eta.prop){
    dist = 1/cophenetic(eb.phylo(tree, tree.ht, eta.prop))
    diag(dist)<-0
    pdist.new = c(dist[i,]%*%Z)
    no0 = which(pd0)
    if(length(no0)){
        if(length(no0)==sum(Z[i,])) likeli = -Inf else 
        likeli = sum(((log(pdist.new)- log(pdist.old))*Z[i,])[-no0])-
            sum((ywU*(pdist.new - pdist.old))[-no0])
    }
    likeli
})

plot(grid1, (aa))





#### Testing all transformations
load('comGMPD.RData')
library(ape)
library(geiger)
source('library.R')
tree <- read.tree('../Data/mammals.tre')
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(com)])

dd<- cophenetic(rescale(tree, "EB", 0))
host.order <- sapply(rownames(dd), function(r) which(r==rownames(com)))
com = com[host.order,]
com =lof(com)

Z= 1*(com>0)


models = list(list(trans = 'BM', grid = seq(0.1,40,length.out= 20)),
    list(trans = 'delta', grid = seq(0.1,4,length.out = 20)),
    list(trans = 'EB', grid = seq(-0.05,0.05,length.out = 20)),
    list(trans = 'kappa', grid = seq(0.1,4,length.out = 20)),
    list(trans = 'lambda', grid = seq(0,1,length.out = 20)),
    list(trans = 'OU', grid = seq(0.1,4,length.out = 20)))

aux = lapply(models, function(r){
    grid = r$grid
    trans = r$trans
    aux =sapply(r$grid, function(eta){
        phy_dist<- cophenetic(rescale(tree, trans, eta))
        dd =1/phy_dist
        diag(dd)<-0
        pdist = dd %*%Z
        P = 1-exp(-pdist)
        roc = rocCurves(Z=Z, Z_cross= Z, P=P, plot=FALSE, bins=400, all=TRUE)
        tb  = ana.table(Z, Z, roc=roc, plot=FALSE)
        cbind(eta=eta, tb=tb)
    })
})


pdf('AUC-parameter-per-Tree-model.pdf')
par(mfrow=c(2,3))
for(i in 1:length(aux)){
    plot(unlist(aux[[i]]['eta',]), unlist(aux[[i]]['tb.auc',]),
         ylab = 'AUC', xlab='parameter',
         type='b' , main=models[[i]]$trans)
}





##### GMPD naive estimate
> aux[,which.max(unlist(aux['tb.auc',]))]
$eta
[1] -0.025

$tb.auc
[1] 84.09

$tb.thresh
[1] 0.716792

$tb.tot.inter
[1] 3966

$tb.hold.out
[1] 0

$tb.pred
[1] NaN

$tb.pred.all
[1] 0.7849218



> aux[,which.max(unlist(aux['tb.auc',]))]
$eta
[1] 0.018

$tb.auc
[1] 78.09

$tb.thresh
[1] 0.005012531

$tb.tot.inter
[1] 4775

$tb.hold.out
[1] 0

$tb.pred
[1] NaN

$tb.pred.all
[1] 0.6425131


com_pa = 1*(com>0)
ZtZ= com_pa%*%t(com_pa)
dist = cophenetic(rescale(tree, 'EB', 0))

pdf('shared-parasites-versus-phy-original.pdf')
plot(dist[upper.tri(dist)], ZtZ[upper.tri(ZtZ)],
     ylab = 'shared parasites',
     xlab = 'pairwise evolutionary distance',
     ylim = c(0,800))
dev.off()



roc = rocCurves(Z=com_pa, Z_cross=com_pa, P=P, all=TRUE, plot=FALSE, bins=400)
ZCV = 1*(roc$P> roc$thresh)
ZtZ= ZCV %*%t(ZCV)
dist = cophenetic(rescale(tree, 'EB', 0))


ZG = 1*(roc.all$P > roc.all$thresh)
ZtZ= ZG %*%t(ZG)
dist = cophenetic(rescale(tree, 'EB', 0))

pdf('shared-parasites-versus-phy-uncertain-noG.pdf')
plot(dist[upper.tri(dist)], ZtZ[upper.tri(ZtZ)],
     ylab = 'shared parasites',
     xlab = 'pairwise evolutionary distance',
     ylim = c(0, 800))
dev.off()
