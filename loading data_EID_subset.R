#!~/bin/R
obj.names = ls()
NOHUMAN = TRUE
TRAIT = FALSE
GMPD_EID_FILTER = FALSE
library(reshape2)
###################
### Loading host-parasite list
# Data from EID only - to be used in the 1st Math paper. 
# Host names standardized to Wilson Reeder 2005 (mammal super tree)
dt <- read.csv("../Data/hp_list_EID_online.csv", as.is=TRUE)

## Cleaning parasite names
## Some parasites have entries with first letter capitalized
dt$parasite <- tolower(dt$parasite)

## This creates 3 host-parasite combinations that are duplicated
## However, counts are very close, so I will remove one of the sets
# dt[which(duplicated(paste0(dt$host,dt$parasite))),]
# dt[which(duplicated(paste0(dt$host,dt$parasite), fromLast=TRUE)),]
dt <- dt[-which(duplicated(paste0(dt$host,dt$parasite))),]

## Keeping only Parasite and Host columns (removing source information)
## If only Publication are used the com has small dimension (43, 260) New   58 436
## If only Sequenences (host, parasite) = ( 530 667) New  410 622
## If Pub + Seq (host, parasite) = (530, 721) New 413 780
dt <- cbind(dt[,c("parasite","host")], count = dt$Sequences + dt$Publications)
# removing "Homo_sapiens" as they have many single-host parasites 
if (NOHUMAN)
    dt <- dt[dt$host!="Homo_sapiens",]

if (GMPD_EID_FILTER){
    pan<-read.csv("../Data/pantheria_WR2005_v11_fMax.csv", header=T)#from Jonathan (BIOL 645 Course)
	dt <- dt[dt$host%in%	pan$bionomial[pan$Order%in%c("Artiodactyla", "Perissodactyla", "Carnivora")],]
	}
###################

### Generating host trait distance matrix
if (TRAIT){
    library(FD)  ##needed when doing full analysis, comparing traits to phylogeny.
    pan<-read.csv("../Data/pantheria_WR2005_v11_fMax.csv", header=T)#from Jonathan (BIOL 645 Course)
    pan <- pan[pan$bionomial %in% dt$host,]
    pan <- pan[,-grep("EXT", colnames(pan))]
    trt <- c(10:38,41)
    pan_trt <- pan[,trt]

                                        # Terrestriality should be converted to 0/1
                                        # TrophicLevel & HabitatBreadth & DietBreadth should be converted to ordered factor
    pan_trt$Terrestriality[pan_trt$Terrestriality==1] <- 0
    pan_trt$Terrestriality[pan_trt$Terrestriality==2] <- 1
    pan_trt$TrophicLevel <- as.ordered(pan_trt$TrophicLevel)
    pan_trt$HabitatBreadth <- as.ordered(pan_trt$HabitatBreadth)
    pan_trt$DietBreadth <- as.ordered(pan_trt$DietBreadth)

    trt_dist <- as.matrix(gowdis(pan_trt))
    ## changing names
    rownames(trt_dist) <- pan$bionomial
    colnames(trt_dist) <- rownames(trt_dist)
# Did this and resulting red/green plot looked bad - almost immediately lost all true 1's
# # 1/distance
# trt_dist= 1/trt_dist
# trt_dist[trt_dist==Inf]<-1

# Instead to make it a similarity matrix, I will take 1-gowdis
# range(trt_dist) # 0 - 0.479
trt_dist <- 1-trt_dist
# range(trt_dist) # 0.52 - 1.00
# trt_dist <- dist_ordering(trt_dist, com)
}
###################


# Subsetting using PanTHERIA to identify clades
pan<-read.csv("../Data/pantheria_WR2005_v11_fMax.csv", header=T)#from Molly
pan <- pan[pan$bionomial %in% dt$host,]
pan <- pan[,-grep("EXT", colnames(pan))]


# Subsetting to only Carnivores => dim(com): 111 x 546
# dt <- dt[dt$host %in% pan$bionomial[pan$Order=="Carnivora"],]

# Subsetting to only Ungulates => dim(com): 97 x 991
# dt <- dt[dt$host %in% pan$bionomial[pan$Order=="Perissodactyla"|pan$Order=="Artiodactyla"],]

# Subsetting to only Primates => dim(com): 114 x 159
# dt <- dt[dt$host %in% pan$bionomial[pan$Order=="Primates"],]

# Subsetting to only Rodents => dim(com): 236 x 442
#dt <- dt[dt$host %in% pan$bionomial[pan$Order=="Rodentia"],]

# Subsetting to only Bats => dim(com): 116 x 60
# dt <- dt[dt$host %in% pan$bionomial[pan$Order=="Chiroptera"],]

## Matrix form where each cell is the FirstCitationDate
## ## Calculating first citation year
## a = regexpr('[0-9]{2,4}+', dt$Citation, perl=TRUE)
## b = substr(dt$Citation, a, a + 3)
## b= as.numeric(b)
## dt$FirstCitationDate<-b
## dt =  dt[c(1,2,3,length(dt), 5:length(dt)-1)]
## #write.csv(dt,file='../Data/GMPD_Feb18_2015_DateFirstCite.csv', row.names = FALSE)

## hp<- unique(data.frame( dt[, c('hostnames', 'parasite')], FirstCitationDate = dt$FirstCitationDate) )
## aux = paste0(hp$hostnames, hp$parasite)
## hp<- unique(hp[, c("hostnames", "parasite")])
## Cite = sapply(paste0(hp$hostnames, hp$parasite), function(r){
## 	min(dt$FirstCitationDate[which(aux==r)])
## 	})
## hp<-cbind(hp, CiteCount = Cite)


### Creating interaction / community matrix
com <- dcast(dt,host ~ parasite, value.var="count", fill=0)
rownames(com)<-com$host
com<-subset(com, select = -host)
com<-as.matrix(com)
# rownames(com)[which(rowSums(com)==max(rowSums(com)))]
dim(com)

## Removing all unnecessary variables
rm(list=grep('(com|pan)', setdiff(ls(), obj.names), invert=TRUE, value=T))

