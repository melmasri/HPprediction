#!~/bin/R
obj.names = ls()
library(reshape2)
### Formating the dataset
#dt <-read.csv('../Data/GMPD_incomplete_Jan23_20155.csv', as.is =TRUE)
dt <-read.csv('../Data/GMPD_Feb18_2015.csv', as.is =TRUE)

# Remove entries with prevalence = 0 : because parasite was tested and not found.
#dt <- subset(dt, Prevalence!=0)
dt <- dt[which(dt$Prevalence!=0),]   # Tested Faster
# Removing NA parasites, this fixed in the final dataset
dt <- dt[!is.na(dt$parasite),]

# Remove parasites recorded to Genus only
dt <- dt[grep("_$",dt$parasite, invert=TRUE),]
dt <- dt[grep("_NA",dt$parasite, invert=TRUE),]
dt <- dt[grep("_sp[.]",dt$parasite, invert=TRUE),]
dt <- dt[grep(" sp[.]",dt$parasite, invert=TRUE),]
dt <- dt[grep("_car ",dt$parasite, invert=TRUE),]

# Removing the following sampling types
#c('Unknown', 'Serology', '')  #ME-MF 20150206: to be removed later(sensitivity).
#table(dt$SamplingType)

###################
# Back to community matrix

# Reformatting to Interaction matrix
# Removing spaces and other non-unique characters
Citation = gsub('(et al|et al.| |-)', '', dt$Citation)

## Calculating unique citation count
hp<- unique(data.frame( dt[, c('hostnames', 'parasite')], Citation=Citation) )
aux = table(paste0(hp$hostnames, hp$parasite))
hp<- unique(hp [, c("hostnames", "parasite")])
CiteCount = sapply(paste0(hp$hostnames, hp$parasite), function(r){
	aux[which(names(aux)==r)]
	})
hp<-cbind(hp, CiteCount)


# Creating interaction / community matrix
com <- dcast(hp,hostnames ~ parasite, value.var='CiteCount', fill=0)
rownames(com)<-com$hostnames
com<-subset(com, select = -hostnames)
com<-as.matrix(com)

## Removing all unnecessary variables
rm(list=grep('com', setdiff(ls(), obj.names), invert=TRUE, value=T))

