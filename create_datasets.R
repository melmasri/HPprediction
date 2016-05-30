## Replicating numbers in Paper
## GMP
rm(list=ls())
source('library.R')
source('loading data_GMPD.R')
source('load_phyDistance.R')

dim(com)
dim(phy_dist)

## Removing
## rows with less than 2 interactions
min.no= 2

## aux = which(rowSums(1*(com>0))<min.no)

## if(length(aux)>0){
##     com = com[-aux, ]
##     phy_dist = phy_dist[-aux,]
##     phy_dist = phy_dist[,-aux]
## }

dim(com)
#empty columns
aux = which(colSums(1*(com>0))<min.no)
if(length(aux)>0)
    com = com[, -aux]
dim(com)

aux = which(rowSums(1*(com>0))==0);aux  #should be none
if(length(aux)>0){
    com = com[-aux,]
    phy_dist= phy_dist[-aux, ]
    phy_dist = phy_dist[, -aux]
}

aux = which(rowSums(1*(com>0))==0);aux  #should be none
aux = which(colSums(1*(com>0))==0);aux  #should be none
aux = which(colSums(1*(com>0))==1);aux  #should be none



com=lof(com)

dim(com)

range( phy_dist)/(max(phy_dist)+1e-2)
phy_dist= phy_dist/(max(phy_dist)+1e-2)
plot_Z(1*(com>0))


save(com, phy_dist, file='comGMPD.RData')

## EID
rm(list=ls())
source('library.R')
source('loading data_EID.R')
source('load_phyDistance.R')
dim(com)
dim(phy_dist)

## Removing
## rows with less than 2 interactions
min.no= 2

dim(com)
#empty columns
aux = which(colSums(1*(com>0))<min.no)
if(length(aux)>0)
    com = com[, -aux]
dim(com)

aux = which(rowSums(1*(com>0))==0);aux  #should be none
if(length(aux)>0){
    com = com[-aux,]
    phy_dist= phy_dist[-aux, ]
    phy_dist = phy_dist[, -aux]
}



aux = which(rowSums(1*(com>0))==0);aux  #should be none
aux = which(colSums(1*(com>0))==0);aux  #should be none
aux = which(colSums(1*(com>0))==1);aux  #should be none

com=lof(com)

dim(com)
dim(phy_dist)
range( phy_dist)/(max(phy_dist)+1e-2)
phy_dist= phy_dist/(max(phy_dist)+1e-2)
plot_Z(1*(com>0))

save(com, phy_dist, file='comEID-PS.RData')


## GMP - subset
rm(list=ls())
source('library.R')
source('loading data_GMPD_subset.R')
source('load_phyDistance.R')
min.no=1
dim(com)
dim(phy_dist)

## Removing
## rows with less than min.no interactions
aux = which(rowSums(1*(com>0))<min.no)
if(length(aux)>0){
    com = com[-aux, ]
    phy_dist = phy_dist[-aux,]
    phy_dist = phy_dist[,-aux]
}
dim(com)
#empty columns
aux = which(colSums(1*(com>0))<min.no)
if(length(aux)>0)
    com = com[, -aux]
dim(com)

aux = which(rowSums(1*(com>0))==0);aux  #should be none
aux = which(colSums(1*(com>0))==0);aux #should be none
com=lof(com)

range( phy_dist)/(max(phy_dist)+1e-2)
phy_dist= phy_dist/(max(phy_dist)+1e-2)
plot_Z(1*(com>0))

save(com, phy_dist,pan, file='comGMPD-year.single.RData')


## EID - subset
rm(list=ls())
source('library.R')
source('loading data_EID_subset.R')
source('load_phyDistance.R')
dim(com)
dim(phy_dist)
min.no=1
## Removing
## rows with less than min.no interactions
aux = which(rowSums(1*(com>0))<min.no)
if(length(aux)>0){
    com = com[-aux, ]
    phy_dist = phy_dist[-aux,]
    phy_dist = phy_dist[,-aux]
}
dim(com)

##empty columns
aux = which(colSums(1*(com>0))<min.no)
if(length(aux)>0){
    com = com[, -aux]
}

dim(com)
aux = which(rowSums(1*(com>0))==0);aux  #should be none
if(length(aux)>0){
    com=com[-aux,]
    phy_dist = phy_dist[-aux,]
    phy_dist = phy_dist[,-aux]
}
aux = which(colSums(1*(com>0))==0);aux #should be none

aux = which(rowSums(1*(com>0))==0);aux  #should be none

com=lof(com)

range( phy_dist)/(max(phy_dist)+1e-2)
phy_dist= phy_dist/(max(phy_dist)+1e-2)
plot_Z(1*(com>0))

save(com, phy_dist, pan, file='comEID-PS-subset.single.RData')
