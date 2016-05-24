#!~/bin/R
#################################3
## loading Distance
## Phylogenetic Distance Matrix
obj.names = ls()
library(ape)
library(geiger)

# Hostnames
if(length(grep('com', ls()))==0)
    stop('No interaction matrix found.')

tree <- read.tree('../Data/mammals.tre')
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(com)])

# Testing all names in hosts in com exist in tree
if(! all(sapply(rownames(com), function(r) r %in% tree$tip.label))) {
		print('Warning! Not all hosts in com exist in tree. Hosts not found in tree will be removed.')
		com <- com[rownames(com)%in%tree$tip.label,]
			}

# testing for lambda
if(length(grep('lambda_phy', ls()))==0){
    print('No lambda is found, default lambda=1')
    lambda_phy=1
}
# Scaling the tree.
phy_dist<- cophenetic(rescale(tree, "lambda", lambda_phy))
# Making sure that no 0 distances off-diagonal
if(sum(diag(phy_dist)==0)!= sum(phy_dist==0))
    stop('Zeros found in off-diagonal of distance matrix.')

## 1/distance

phy_dist= 1/phy_dist
phy_dist[phy_dist==Inf]<-0

#Ordering the distance matrix with the corresponding interaction matrix
phy_dist = dist_ordering(phy_dist, com)


rm(list=grep('phy_dist', setdiff(ls(), obj.names), invert=TRUE, value=T))
## Testiing
## lambda_phy =0.8
## range(phy_dist[upper.tri(phy_dist)])
## d= max(1/phy_dist);d
## range(phy_dist)
## summary(range(phy_dist[upper.tri(phy_dist)])  )
## max(phy_dist)


