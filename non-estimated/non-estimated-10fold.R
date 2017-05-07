############################################
## Script to run 10fold cross validation
## Non-estimated models

rm(list=ls())

load("comGMPD.RData")

# Convert count to binary
com[com>1] <- 1
com <- unname(com)

#######################
## set the correct prior.
## Summary statistics
cnames = colnames(com)
rnames = rownames(com)
com = unname(com)
phy_dist = unname(phy_dist)

source('../library.R', local=TRUE)
pairs = cross.validate.fold(1*(com>0), n=5)
tot.gr = length(unique(pairs[,'gr']))

est <- function(x, pairs, Z, dist_type, phy_dist){
        source('../library.R', local=TRUE)
        source('../gen.R', local=TRUE)

        com_paCross = Z
        com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
            
        if(dist_type == 'phy'){
            dist = phy_dist
            P = 1-exp(-(dist)%*%com_paCross) 
        }

        if(dist_type == 'ztz'){

            dist = com_paCross%*%t(com_paCross)
            diag(dist) = 0
            dist = dist/max(dist)
            dist = 1-(1/dist)
            P = 1-exp(-(dist)%*%com_paCross) 
        }

        if(dist_type == 'uniform'){

            dist = matrix(1, nrow=nrow(com), ncol=nrow(com))
            diag(dist) = 0
            P = 1-exp(-(dist)%*%com_paCross) 
        }

        if(dist_type == 'ztz_phy'){
            phy = phy_dist
            dist = com_paCross%*%t(com_paCross)
            diag(dist) = 0
            dist = dist/max(dist)
            dist = 1-(1/dist)
            P = 1-exp(-((dist)%*%com_paCross) + (phy)%*%com_paCross)) 
        }

        roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
        tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
        roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
        tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
        list(tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
        }

res_unif = lapply(1:tot.gr, est, pairs=pairs, Z = com, dist_type = "uniform", phy_dist=phy_dist) 
res_ztz = lapply(1:tot.gr, est, pairs=pairs, Z = com, dist_type = "ztz", phy_dist=phy_dist) 
res_phy = lapply(1:tot.gr, est, pairs=pairs, Z = com, dist_type = "phy", phy_dist=phy_dist) 


# Averaging results over each fold
m.auc = lapply(list(res_unif=res_unif,res_ztz=res_ztz,res_phy=res_phy), function(x) mean(sapply(x, function(r) r$tb$auc)))
m.thresh = lapply(list(res_unif=res_unif,res_ztz=res_ztz,res_phy=res_phy), function(x) mean(sapply(x, function(r) r$tb$thres)))
m.pred = lapply(list(res_unif=res_unif,res_ztz=res_ztz,res_phy=res_phy), function(x) mean(sapply(x, function(r) r$tb$pred)))
m.hold.out = lapply(list(res_unif=res_unif,res_ztz=res_ztz,res_phy=res_phy), function(x) mean(sapply(x, function(r) r$tb$hold.out)))

m.auc.all = lapply(list(res_unif=res_unif,res_ztz=res_ztz,res_phy=res_phy), function(x) mean(sapply(x, function(r) r$tb.all$auc)))
m.thresh.all = lapply(list(res_unif=res_unif,res_ztz=res_ztz,res_phy=res_phy), function(x) mean(sapply(x, function(r) r$tb.all$thres)))
m.pred.all = lapply(list(res_unif=res_unif,res_ztz=res_ztz,res_phy=res_phy), function(x) mean(sapply(x, function(r) r$tb.all$pred)))
m.hold.out.all = lapply(list(res_unif=res_unif,res_ztz=res_ztz,res_phy=res_phy), function(x) mean(sapply(x, function(r) r$tb.all$hold.out)))


# ztz is keeping more zeros
# phylogeny is pushign some additional zeros to ones

# if(SAVE_PARAM)
#     save.image(file = 'param.RData')




##################################################
