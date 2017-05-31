# tree transformations

require(ape)
# require(geiger)

set.seed(20150508)

# Herc's transformation
tree<-rcoal(250)
plot(tree)
dist <- cophenetic(tree)


# Pagel's Lambda
lambdas <- seq(0,1,0.05)
lambda_trees <- list(NULL)

for (i in seq_along(lambdas)){
	lambda_trees[i] <- list(rescale(tree, "lambda", lambdas[i]))
}


par(mfrow=c(1,4))
plot(lambda_trees[[1]])
title("lambda: 0.0")
plot(lambda_trees[[7]])
title("lambda: 0.3")
plot(lambda_trees[[15]])
title("lambda: 0.7")
plot(lambda_trees[[21]])
title("lambda: 1")


# Pagel's Delta
deltas <- seq(0,4,0.5)
delta_trees <- list(NULL)

for (i in seq_along(deltas)){
	delta_trees[i] <- list(rescale(tree, "delta", deltas[i]))
}

par(mfrow=c(1,4))
plot(delta_trees[[1]])
title("delta: 0.0")
plot(delta_trees[[2]])
title("delta: 0.5")
plot(delta_trees[[5]])
title("delta: 2")
plot(delta_trees[[9]])
title("delta: 4")


# Pagel's kappa
kappas <- seq(0,3,0.5)
kappa_trees <- list(NULL)

for (i in seq_along(kappas)){
	kappa_trees[i] <- list(rescale(tree, "kappa", kappas[i]))
}

par(mfrow=c(1,4))
plot(kappa_trees[[1]])
title("kappa: 0.0")
plot(kappa_trees[[2]])
title("kappa: 0.5")
plot(kappa_trees[[3]])
title("kappa: 1.0")
plot(kappa_trees[[5]])
title("kappa: 2.0")


# OU model
alphas <- seq(0,20,0.5)
alpha_trees <- list(NULL)

for (i in seq_along(alphas)){
	alpha_trees[i] <- list(rescale(tree, "OU", alphas[i]))
}

par(mfrow=c(1,4))
plot(alpha_trees[[2]])
title("alpha: 0.5")
plot(alpha_trees[[5]])
title("alpha: 2.0")
plot(alpha_trees[[21]])
title("alpha: 10.0")
plot(alpha_trees[[41]])
title("alpha: 20")


# EB model
# Rate (a) of zero returns original tree
rates <- seq(-1,2,0.1)
rate_trees <- list(NULL)

for (i in seq_along(rates)){
	rate_trees[i] <- list(rescale(tree, "EB", rates[i]))
}

par(mfrow=c(1,4))
plot(rate_trees[[1]])
title("rate: -1.0")
plot(rate_trees[[11]])
title("rate: 0.0")
plot(rate_trees[[21]])
title("rate: 1.0")
plot(rate_trees[[31]])
title("rate: 2.0")

par(mfrow=c(1,1))
plot(rate_trees[[31]])

is.ultrametric(rate_trees[[31]])
is.ultrametric(rate_trees[[1]])



############################################
############################################
## Script to run 5 fold cross validation
## Non-estimated models 
## with varying tree transformations

# GMPD
load("comGMPD.RData")

# Convert count to binary
com[com>1] <- 1

#!~/bin/R
#################################3
## loading Distance
## Phylogenetic Distance Matrix
obj.names = ls()
library(ape)
library(geiger)

tree <- read.tree('../../Data/mammals.tre')
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(com)])

# com <- unname(com)

#######################
## set the correct prior.
## Summary statistics
cnames = colnames(com)
rnames = rownames(com)
# com = unname(com)
# phy_dist = unname(phy_dist)

source('../library.R', local=TRUE)
pairs = cross.validate.fold(1*(com>0), n=5)
tot.gr = length(unique(pairs[,'gr']))

est_transform <- function(x, pairs, Z, tree, trans_type, trans_param){
        # source('../library.R', local=TRUE)
        # source('../gen.R', local=TRUE)

        com_paCross = Z
        com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
        

		# Scaling the tree.
		phy_dist<- cophenetic(rescale(tree, trans_type, trans_param))

		## 1/distance
		phy_dist= 1/phy_dist
		phy_dist[phy_dist==Inf]<-0

		#Ordering the distance matrix with the corresponding interaction matrix
		phy_dist = dist_ordering(phy_dist, com)

        P = 1-exp(-(phy_dist)%*%com_paCross) 

        roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
        tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
        roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
        tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
        list(tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
        }

par(mfrow=c(1,1))

# Pagel's Lambda
lambdas <- seq(0,1,0.05)
res_lambda <- list(NULL)
for (i in seq_along(lambdas)){
    res_lambda[i] <- list(lapply(1:tot.gr, est_transform, pairs=pairs, Z = com, tree=tree, trans_type="lambda", trans_param=lambdas[i])) 
}

# Pagel's Delta 
deltas <- seq(0.1,40,1.5)
res_delta <- list(NULL)
for (i in seq_along(deltas)){
    res_delta[i] <- list(lapply(1:tot.gr, est_transform, pairs=pairs, Z = com, tree=tree, trans_type="delta", trans_param=deltas[i])) 
}

# Pagel's Kappa
kappas <- c(seq(0,1.8,0.2),seq(2,6,0.4))
length(kappas)
res_kappa <- list(NULL)
for (i in seq_along(kappas)){
    res_kappa[i] <- list(lapply(1:tot.gr, est_transform, pairs=pairs, Z = com, tree=tree, trans_type="kappa", trans_param=kappas[i])) 
}


# OU model
alphas <- seq(0.1,10,0.2)
length(alphas)
res_ou <- list(NULL)
for (i in seq_along(alphas)){
    res_ou[i] <- list(lapply(1:tot.gr, est_transform, pairs=pairs, Z = com, tree=tree, trans_type="OU", trans_param=alphas[i])) 
}


# EB model
# Rate (a) of zero returns original tree
rates <- seq(-0.1,0.05,0.005)
length(rates)
res_eb <- list(NULL)
for (i in seq_along(rates)){
    res_eb[i] <- list(lapply(1:tot.gr, est_transform, pairs=pairs, Z = com, tree=tree, trans_type="EB", trans_param=rates[i])) 
}

ls()
# save.image("GMPD_tree_transforms.RData")
load("GMPD_tree_transforms.RData")

# Results table
# res <- list(lambda=res_lambda, delta=res_delta, kappa=res_kappa, ou=res_ou, eb=res_eb)
res <- c(res_lambda, res_delta, res_kappa, res_ou, res_eb)
trans_name <- c("lambda", "delta", "kappa", "ou", "eb")

df_trans_gmpd <- NULL

for (i in seq_along(res)){

	# name <- trans_name[i]

	m.auc = lapply(res[i], function(x) mean(sapply(x, function(r) r$tb$auc)))
	m.thresh = lapply(res[i], function(x) mean(sapply(x, function(r) r$tb$thres)))
	m.pred = lapply(res[i], function(x) mean(sapply(x, function(r) r$tb$pred)))
	m.hold.out = lapply(res[i], function(x) mean(sapply(x, function(r) r$tb$hold.out)))

	m.auc.all = lapply(res[i], function(x) mean(sapply(x, function(r) r$tb.all$auc)))
	m.thresh.all = lapply(res[i], function(x) mean(sapply(x, function(r) r$tb.all$thres)))
	m.pred.all = lapply(res[i], function(x) mean(sapply(x, function(r) r$tb.all$pred)))
	m.hold.out.all = lapply(res[i], function(x) mean(sapply(x, function(r) r$tb.all$hold.out)))
	
	df <- data.frame(Database="GMPD", m.auc=unlist(m.auc),m.thresh=unlist(m.thresh),m.pred=unlist(m.pred),m.hold.out=unlist(m.hold.out),
                 m.auc.all=unlist(m.auc.all),m.thresh.all=unlist(m.thresh.all),m.pred.all=unlist(m.pred.all),m.hold.out.all=unlist(m.hold.out.all))
	df_trans_gmpd <- rbind(df_trans_gmpd, df)
}

df_trans_gmpd <- cbind(value = c(lambdas,deltas,kappas,alphas,rates),df_trans_gmpd)
df_trans_gmpd <- cbind(trans = c(rep("lambda", length(lambdas)),
								 rep("delta", length(deltas)),
								 rep("kappa", length(kappas)),
								 rep("ou", length(alphas)),
								 rep("eb", length(rates))),df_trans_gmpd)

str(df_trans_gmpd)

# Plotting
for (i in 1:length(unique(df_trans_gmpd$trans))) {
	transf <- unique(df_trans_gmpd$trans)[i]
	pdf(paste("GMPD",transf,"5foldCV.pdf", sep="_"))
	dat <- df_trans_gmpd[df_trans_gmpd$trans==transf,]
	plot(dat$m.auc ~ dat$value, xlab=paste(transf), ylab="Mean AUC across 5 folds", main=paste0("AUC for phylogenetic distance transformed by ",transf))
	lines(dat$m.auc ~ dat$value)
	with(dat, abline(v=value[m.auc==max(m.auc)], col=2, lty=2))
	dev.off()
}



# # Averaging results over each fold
# m.auc = lapply(res_lambda, function(x) mean(sapply(x, function(r) r$tb$auc)))
# plot(unlist(m.auc)~lambdas)

# m.auc = lapply(res_delta, function(x) mean(sapply(x, function(r) r$tb$auc)))
# plot(unlist(m.auc)~deltas)

# m.auc = lapply(res_kappa, function(x) mean(sapply(x, function(r) r$tb$auc)))
# plot(unlist(m.auc)~kappas)

# m.auc = lapply(res_ou, function(x) mean(sapply(x, function(r) r$tb$auc)))
# plot(unlist(m.auc)~alphas)

# m.auc = lapply(res_eb, function(x) mean(sapply(x, function(r) r$tb$auc)))
# plot(unlist(m.auc)~rates)



### EID

load("comEID-PS.RData")


# Convert count to binary
com[com>1] <- 1

#!~/bin/R
#################################3
## loading Distance
## Phylogenetic Distance Matrix
obj.names = ls()
library(ape)
library(geiger)

tree <- read.tree('../../Data/mammals.tre')
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(com)])

# com <- unname(com)

#######################
## set the correct prior.
## Summary statistics
cnames = colnames(com)
rnames = rownames(com)
# com = unname(com)
# phy_dist = unname(phy_dist)

source('../library.R', local=TRUE)
pairs = cross.validate.fold(1*(com>0), n=5)
tot.gr = length(unique(pairs[,'gr']))

# par(mfrow=c(1,1))

# Pagel's Lambda
lambdas <- seq(0,1,0.05)
res_lambda <- list(NULL)
for (i in seq_along(lambdas)){
    res_lambda[i] <- list(lapply(1:tot.gr, est_transform, pairs=pairs, Z = com, tree=tree, trans_type="lambda", trans_param=lambdas[i])) 
}

# Pagel's Delta 
deltas <- seq(0.1,40,1.5)
res_delta <- list(NULL)
for (i in seq_along(deltas)){
    res_delta[i] <- list(lapply(1:tot.gr, est_transform, pairs=pairs, Z = com, tree=tree, trans_type="delta", trans_param=deltas[i])) 
}

# Pagel's Kappa
kappas <- c(seq(0,1.8,0.2),seq(2,6,0.4))
length(kappas)
res_kappa <- list(NULL)
for (i in seq_along(kappas)){
    res_kappa[i] <- list(lapply(1:tot.gr, est_transform, pairs=pairs, Z = com, tree=tree, trans_type="kappa", trans_param=kappas[i])) 
}


# OU model
alphas <- seq(0.1,10,0.2)
length(alphas)
res_ou <- list(NULL)
for (i in seq_along(alphas)){
    res_ou[i] <- list(lapply(1:tot.gr, est_transform, pairs=pairs, Z = com, tree=tree, trans_type="OU", trans_param=alphas[i])) 
}


# EB model
# Rate (a) of zero returns original tree
rates <- seq(-0.1,0.05,0.005)
length(rates)
res_eb <- list(NULL)
for (i in seq_along(rates)){
    res_eb[i] <- list(lapply(1:tot.gr, est_transform, pairs=pairs, Z = com, tree=tree, trans_type="EB", trans_param=rates[i])) 
}

ls()
# save.image("EID_tree_transforms.RData")
# load("EID_tree_transforms.RData")

# Results table
# res <- list(lambda=res_lambda, delta=res_delta, kappa=res_kappa, ou=res_ou, eb=res_eb)
res <- c(res_lambda, res_delta, res_kappa, res_ou, res_eb)
trans_name <- c("lambda", "delta", "kappa", "ou", "eb")

df_trans_eid <- NULL

for (i in seq_along(res)){

	# name <- trans_name[i]

	m.auc = lapply(res[i], function(x) mean(sapply(x, function(r) r$tb$auc)))
	m.thresh = lapply(res[i], function(x) mean(sapply(x, function(r) r$tb$thres)))
	m.pred = lapply(res[i], function(x) mean(sapply(x, function(r) r$tb$pred)))
	m.hold.out = lapply(res[i], function(x) mean(sapply(x, function(r) r$tb$hold.out)))

	m.auc.all = lapply(res[i], function(x) mean(sapply(x, function(r) r$tb.all$auc)))
	m.thresh.all = lapply(res[i], function(x) mean(sapply(x, function(r) r$tb.all$thres)))
	m.pred.all = lapply(res[i], function(x) mean(sapply(x, function(r) r$tb.all$pred)))
	m.hold.out.all = lapply(res[i], function(x) mean(sapply(x, function(r) r$tb.all$hold.out)))
	
	df <- data.frame(Database="EID", m.auc=unlist(m.auc),m.thresh=unlist(m.thresh),m.pred=unlist(m.pred),m.hold.out=unlist(m.hold.out),
                 m.auc.all=unlist(m.auc.all),m.thresh.all=unlist(m.thresh.all),m.pred.all=unlist(m.pred.all),m.hold.out.all=unlist(m.hold.out.all))
	df_trans_eid <- rbind(df_trans_eid, df)
}

df_trans_eid <- cbind(value = c(lambdas,deltas,kappas,alphas,rates),df_trans_eid)
df_trans_eid <- cbind(trans = c(rep("lambda", length(lambdas)),
								 rep("delta", length(deltas)),
								 rep("kappa", length(kappas)),
								 rep("ou", length(alphas)),
								 rep("eb", length(rates))),df_trans_eid)

str(df_trans_eid)

# Plotting
for (i in 1:length(unique(df_trans_eid$trans))) {
	transf <- unique(df_trans_eid$trans)[i]
	pdf(paste("EID",transf,"5foldCV.pdf", sep="_"))
	dat <- df_trans_eid[df_trans_eid$trans==transf,]
	plot(dat$m.auc ~ dat$value, xlab=paste(transf), ylab="Mean AUC across 5 folds", main=paste0("AUC for phylogenetic distance transformed by ",transf))
	lines(dat$m.auc ~ dat$value)
	with(dat, abline(v=value[m.auc==max(m.auc)], col=2, lty=2))
	dev.off()
}


### table

# top m.auc
df_trans_eid[df_trans_eid$m.auc==max(df_trans_eid$m.auc),]
df_trans_gmpd[df_trans_gmpd$m.auc==max(df_trans_gmpd$m.auc),]

# top m.auc.all
df_trans_eid[df_trans_eid$m.auc.all==max(df_trans_eid$m.auc.all),]
df_trans_gmpd[df_trans_gmpd$m.auc.all==max(df_trans_gmpd$m.auc.all),]

df_full <- rbind(df_trans_gmpd, df_trans_eid)

library(xtable)
tab <- xtable(df_full)
print(tab, file="tree_trans_table.tex")

# EB model
# Rate (a) of zero returns original tree
rates <- c(-0.035, -0.02, 0, 0.02)
eb_trees <- list(NULL)
for (i in seq_along(rates)){
	eb_trees[i] <- list(rescale(tree, "EB", rates[i]))
}

par(mfrow=c(1,4))
plot(eb_trees[[1]], tip.label=FALSE)
title("rate: -0.035")
plot(eb_trees[[2]], tip.label=FALSE)
title("rate: -0.02")
plot(eb_trees[[3]], tip.label=FALSE)
title("rate: 0.00")
plot(eb_trees[[4]], tip.label=FALSE)
title("rate: 0.02")

# interpretation...
df_full[df_full$m.auc>87,]
df_full[df_full$m.auc.all>82,]
