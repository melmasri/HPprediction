#!~/bin/R
obj.names = ls()
library(reshape2)

### Download and load data
if(!file.exists("example/GMPD_main.csv")){
	require(fulltext)
	# GMPD v2.0
	# Stephens et al. 2016 Ecology
	dt <- read.csv(unzip(ft_get_si("10.1002/ecy.1799", 1, "wiley"), exdir="example/", "GMPD_datafiles/GMPD_main.csv", junkpaths=TRUE), as.is=TRUE)
} else dt <- read.csv("example/GMPD_main.csv", as.is=TRUE)

### Formating the dataset
# Subsetting to Carnivore + Ungulate subsections
dt <- dt[dt$Group %in% c("carnivores","ungulates"),]

# Removing parasites not reported to species
Sys.setlocale('LC_ALL','C') 
dt <- dt[grep("sp[.]",dt$ParasiteCorrectedName, invert=TRUE),]
dt <- dt[grep("ABOLISHED",dt$ParasiteCorrectedName, invert=TRUE),]
dt <- dt[grep("no binomial name",dt$ParasiteCorrectedName, invert=TRUE),]
dt <- dt[grep("not identified to genus",dt$ParasiteCorrectedName, invert=TRUE),]
dt <- dt[grep("SPLITTED in ICTV",dt$ParasiteCorrectedName, invert=TRUE),]
dt <- dt[grep("Diphyllobothrium sp",dt$ParasiteCorrectedName, invert=TRUE),]

# Removing entires with "no binomial name" for hosts
dt <- dt[dt$HostCorrectedName!="no binomial name",]

# Formatting host names to match phylogeny
dt$HostCorrectedName <- gsub(" ","_", dt$HostCorrectedName)

# Remove entries with prevalence = 0 : because ParasiteCorrectedName was tested and not found.
dt <- dt[which(dt$Prevalence!=0),]


### Creating interaction matrix
# Removing spaces and other non-unique characters in Citation
Citation <- gsub('(et al|et al.| |-)', '', dt$Citation)

# Calculating unique citation count
hp <- unique(data.frame( dt[, c('HostCorrectedName', 'ParasiteCorrectedName')], Citation=Citation))
aux <- table(paste0(hp$HostCorrectedName, hp$ParasiteCorrectedName))
hp <- unique(hp [, c("HostCorrectedName", "ParasiteCorrectedName")])
CiteCount <- sapply(paste0(hp$HostCorrectedName, hp$ParasiteCorrectedName), function(r){
	aux[which(names(aux)==r)]
	})
hp <- cbind(hp, CiteCount)

# Creating interaction / community matrix
com <- dcast(hp,HostCorrectedName ~ ParasiteCorrectedName, value.var='CiteCount', fill=0)
rownames(com) <- com$HostCorrectedName
com <- subset(com, select = -HostCorrectedName)
com <- as.matrix(com)

### Removing all unnecessary variables
rm(list=grep('com', setdiff(ls(), obj.names), invert=TRUE, value=T))
