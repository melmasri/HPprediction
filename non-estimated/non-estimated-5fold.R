############################################
## Script to run 10fold cross validation
## Non-estimated models

rm(list=ls())


# GMPD

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
        # source('../library.R', local=TRUE)
        # source('../gen.R', local=TRUE)

        com_paCross = Z
        com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
            
        if(dist_type == 'phy'){
            dist = phy_dist
            P = 1-exp(-(dist)%*%com_paCross) 
        }

        if(dist_type == 'ztz'){

            dist = com_paCross%*%t(com_paCross)
            # dist = 1-(1/dist)
            dist = dist + 1e-4
            diag(dist) = 0
            # dist[is.infinite(dist)] = 0
            dist = dist/ncol(Z)
            P = 1-exp(-(dist)%*%com_paCross) 
        }

        if(dist_type == 'uniform'){

            dist = matrix(1, nrow=nrow(com), ncol=nrow(com))
            diag(dist) = 0
            P = 1-exp(-(dist)%*%com_paCross) 
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

# Results table
df_gmpd <- data.frame(Database="GMPD",m.auc=unlist(m.auc),m.thresh=unlist(m.thresh),m.pred=unlist(m.pred),m.hold.out=unlist(m.hold.out),
                 m.auc.all=unlist(m.auc.all),m.thresh.all=unlist(m.thresh.all),m.pred.all=unlist(m.pred.all),m.hold.out.all=unlist(m.hold.out.all))
df_gmpd
# saveRDS(df_gmpd, "gmpd_non-est_5fold.rds")


# ROC plots

# From paper-code.R

nn.files <- list(res_unif=res_unif,res_ztz=res_ztz,res_phy=res_phy)

test <- nn.files[1]
names(test)
str(test)

gres= lapply(nn.files, function(f){
    # #a =regexpr(paste0('-[A-Za-z0-9]*',data,'-.*h[0-9]{2}'), f, perl=TRUE)
    # a =regexpr(paste0('[A-Za-z-0-9-]*10foldCV-[A-za-z0-9-]*',data), f, perl=TRUE)
    # name = substr(f, a, a + attributes(a)$match.length-1)
    # name =sub('10foldCV-', '', name)
    # print(name)
    # load(f)
    FPR = rowMeans(sapply(f, function(r) r$FPR))
    TPR = rowMeans(sapply(f, function(r) r$TPR))
    ## if(!grepl('-NN[0-9]{2}', f)){       
    m.auc= sapply(f, function(r) r$tb$auc)
    m.thresh=sapply(f, function(r) r$tb$thres)
    m.pred=sapply(f,function(r) r$tb$pred)
    m.hold.out=sapply(f, function(r) r$tb$hold.out)
    ## }else{
    ##     m.auc= sapply(res, function(r) r$auc)
    ##     m.thresh=sapply(res, function(r) r$thres)
    ##     m.pred=sapply(res,function(r) r$pred)
    ##     m.hold.out=sapply(res, function(r) r$hold.out)
    ## }
    # name = paste0(name, round(mean(m.auc),2),'-',round(mean(m.pred),2))
    TB = data.frame(m.auc, m.thresh, m.pred, m.hold.out)
    # list(name=name, graph=cbind(FPR, TPR), ana=TB)
    list(graph=cbind(FPR, TPR), ana=TB)
})
str(gres)

pdf('GMPD-non-est-5foldCV.pdf')
# gnames = names(nn.files)
gnames = c("Uniform","ZtZ","Phylogeny")
# gnames= c('LS-network: full model', 'Nearest-neighbour', 'LS-network: phylogeny-only','LS-network: affinity-only','LS-network: weighted-by-counts')
gcol = c('red', 'blue', 'darkgreen')
glty = c(3,2,4)
#gpch = c('+', '*', 'o', 6, 7)
gpch = c(3, 5, 1)
t = 'b'
glwd=3
i= 1
plot(gres[[i]]$graph, xlab='1-specificity', ylab = 'sensitivity', type =t, col=gcol[i], main = 'ROC Curve', xlim = c(0,1), ylim = c(0,1), lty=glty[i], lwd=glwd, pch=gpch[i])
i=2
lines(gres[[i]]$graph,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
i =3
lines(gres[[i]]$graph,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
legend('bottomright', legend = gnames, col=gcol, lty=glty,lwd=2, pch = gpch)
# legend('bottomright', legend = gnames, col=gcol, lty=glty,lwd=2, pch = gpch)
dev.off()


### TRYING PHYLOGENY TRANSFORMATION

est_eta <- function(x, pairs, Z, dist_type, phy_dist, eta){
        # source('../library.R', local=TRUE)
        # source('../gen.R', local=TRUE)

        com_paCross = Z
        com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
            
        if(dist_type == 'phy'){
            dist = phy_dist
            P = 1-exp(-(dist^eta) %*% com_paCross)
            # P = 1-exp(-(dist)%*%com_paCross) 
        }

        if(dist_type == 'ztz'){

            dist = com_paCross%*%t(com_paCross)
            # dist = 1-(1/dist)
            dist = dist + 1e-4
            diag(dist) = 0
            # dist[is.infinite(dist)] = 0
            dist = dist/ncol(Z)
            P = 1-exp(-(dist)%*%com_paCross) 
        }

        if(dist_type == 'uniform'){

            dist = matrix(1, nrow=nrow(com), ncol=nrow(com))
            diag(dist) = 0
            P = 1-exp(-(dist)%*%com_paCross) 
        }

        roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=FALSE)
        tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
        roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400, all=TRUE)
        tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
        list(tb = tb, tb.all = tb.all, FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR, FPR = roc$roc$FPR, TPR=roc$roc$TPR)
        }


etas <- seq(0,2,0.1)

res_etas <- list(NULL)

for (i in seq_along(etas)){

    res_etas[i] <- list(lapply(1:tot.gr, est_eta, pairs=pairs, Z = com, dist_type = "phy", phy_dist=phy_dist, eta=etas[i])) 
}

length(res_etas)

# Averaging results over each fold
m.auc = lapply(res_etas, function(x) mean(sapply(x, function(r) r$tb$auc)))
m.thresh = lapply(res_etas, function(x) mean(sapply(x, function(r) r$tb$thres)))
m.pred = lapply(res_etas, function(x) mean(sapply(x, function(r) r$tb$pred)))
m.hold.out = lapply(res_etas, function(x) mean(sapply(x, function(r) r$tb$hold.out)))

m.auc.all = lapply(res_etas, function(x) mean(sapply(x, function(r) r$tb.all$auc)))
m.thresh.all = lapply(res_etas, function(x) mean(sapply(x, function(r) r$tb.all$thres)))
m.pred.all = lapply(res_etas, function(x) mean(sapply(x, function(r) r$tb.all$pred)))
m.hold.out.all = lapply(res_etas, function(x) mean(sapply(x, function(r) r$tb.all$hold.out)))

# Results table
df_etas_gmpd <- data.frame(Database="GMPD",eta = etas, m.auc=unlist(m.auc),m.thresh=unlist(m.thresh),m.pred=unlist(m.pred),m.hold.out=unlist(m.hold.out),
                 m.auc.all=unlist(m.auc.all),m.thresh.all=unlist(m.thresh.all),m.pred.all=unlist(m.pred.all),m.hold.out.all=unlist(m.hold.out.all))
rownames(df_etas_gmpd) <- NULL
df_etas_gmpd

# Plotting
pdf("gmpd_phy-dist-eta_5foldCV.pdf")
plot(df_etas_gmpd$m.auc ~ df_etas_gmpd$eta, xlab="eta",ylab="Mean AUC across 5 folds", main="AUC for phylogenetic distance transformed by eta")
lines(df_etas_gmpd$m.auc ~ df_etas_gmpd$eta)
with(df_etas_gmpd, abline(v=eta[m.auc==max(m.auc)], col=2, lty=2))
dev.off()



##################################################


# EID

load("comEID-PS.RData")

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
        # source('../library.R', local=TRUE)
        # source('../gen.R', local=TRUE)

        com_paCross = Z
        com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
            
        if(dist_type == 'phy'){
            dist = phy_dist
            P = 1-exp(-(dist)%*%com_paCross) 
        }

        if(dist_type == 'ztz'){

            dist = com_paCross%*%t(com_paCross)
            # dist = 1-(1/dist)
            dist = dist + 1e-4
            diag(dist) = 0
            # dist[is.infinite(dist)] = 0
            dist = dist/ncol(Z)
            P = 1-exp(-(dist)%*%com_paCross) 
        }

        if(dist_type == 'uniform'){

            dist = matrix(1, nrow=nrow(com), ncol=nrow(com))
            diag(dist) = 0
            P = 1-exp(-(dist)%*%com_paCross) 
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


# Results table
df_eid <- data.frame(Database="EID",m.auc=unlist(m.auc),m.thresh=unlist(m.thresh),m.pred=unlist(m.pred),m.hold.out=unlist(m.hold.out),
                 m.auc.all=unlist(m.auc.all),m.thresh.all=unlist(m.thresh.all),m.pred.all=unlist(m.pred.all),m.hold.out.all=unlist(m.hold.out.all))
df <- cbind(Model=rep(c("Uniform","ZtZ","Phylogeny"),2),rbind(df_gmpd,df_eid))
rownames(df) <- NULL

library(xtable)
tab <- xtable(df)

saveRDS(df,"non-est-5foldCV.rds")
print(tab, file="non-est-5foldCV.tex")


# ROC plots
# From paper-code.R

nn.files <- list(res_unif=res_unif,res_ztz=res_ztz,res_phy=res_phy)

gres= lapply(nn.files, function(f){
    # #a =regexpr(paste0('-[A-Za-z0-9]*',data,'-.*h[0-9]{2}'), f, perl=TRUE)
    # a =regexpr(paste0('[A-Za-z-0-9-]*10foldCV-[A-za-z0-9-]*',data), f, perl=TRUE)
    # name = substr(f, a, a + attributes(a)$match.length-1)
    # name =sub('10foldCV-', '', name)
    # print(name)
    # load(f)
    FPR = rowMeans(sapply(f, function(r) r$FPR))
    TPR = rowMeans(sapply(f, function(r) r$TPR))
    ## if(!grepl('-NN[0-9]{2}', f)){       
    m.auc= sapply(f, function(r) r$tb$auc)
    m.thresh=sapply(f, function(r) r$tb$thres)
    m.pred=sapply(f,function(r) r$tb$pred)
    m.hold.out=sapply(f, function(r) r$tb$hold.out)
    ## }else{
    ##     m.auc= sapply(res, function(r) r$auc)
    ##     m.thresh=sapply(res, function(r) r$thres)
    ##     m.pred=sapply(res,function(r) r$pred)
    ##     m.hold.out=sapply(res, function(r) r$hold.out)
    ## }
    # name = paste0(name, round(mean(m.auc),2),'-',round(mean(m.pred),2))
    TB = data.frame(m.auc, m.thresh, m.pred, m.hold.out)
    # list(name=name, graph=cbind(FPR, TPR), ana=TB)
    list(graph=cbind(FPR, TPR), ana=TB)
})
str(gres)


pdf("EID_non-est-5foldCV.pdf")
# gnames = names(nn.files)
gnames = c("Uniform","ZtZ","Phylogeny")
# gnames= c('LS-network: full model', 'Nearest-neighbour', 'LS-network: phylogeny-only','LS-network: affinity-only','LS-network: weighted-by-counts')
gcol = c('red', 'blue', 'darkgreen')
glty = c(3,2,4)
#gpch = c('+', '*', 'o', 6, 7)
gpch = c(3, 5, 1)
t = 'b'
glwd=3
i= 1
plot(gres[[i]]$graph, xlab='1-specificity', ylab = 'sensitivity', type =t, col=gcol[i], main = 'ROC Curve', xlim = c(0,1), ylim = c(0,1), lty=glty[i], lwd=glwd, pch=gpch[i])
i=2
lines(gres[[i]]$graph,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
i =3
lines(gres[[i]]$graph,type=t, col=gcol[i],lty=glty[i], lwd=glwd,pch=gpch[i])
legend('bottomright', legend = gnames, col=gcol, lty=glty,lwd=2, pch = gpch)
# legend('bottomright', legend = gnames, col=gcol, lty=glty,lwd=2, pch = gpch)
dev.off()


### TRYING PHYLOGENY TRANSFORMATION

etas <- seq(0,2,0.1)

res_etas <- list(NULL)

for (i in seq_along(etas)){

    res_etas[i] <- list(lapply(1:tot.gr, est_eta, pairs=pairs, Z = com, dist_type = "phy", phy_dist=phy_dist, eta=etas[i])) 
}

length(res_etas)

# Averaging results over each fold
m.auc = lapply(res_etas, function(x) mean(sapply(x, function(r) r$tb$auc)))
m.thresh = lapply(res_etas, function(x) mean(sapply(x, function(r) r$tb$thres)))
m.pred = lapply(res_etas, function(x) mean(sapply(x, function(r) r$tb$pred)))
m.hold.out = lapply(res_etas, function(x) mean(sapply(x, function(r) r$tb$hold.out)))

m.auc.all = lapply(res_etas, function(x) mean(sapply(x, function(r) r$tb.all$auc)))
m.thresh.all = lapply(res_etas, function(x) mean(sapply(x, function(r) r$tb.all$thres)))
m.pred.all = lapply(res_etas, function(x) mean(sapply(x, function(r) r$tb.all$pred)))
m.hold.out.all = lapply(res_etas, function(x) mean(sapply(x, function(r) r$tb.all$hold.out)))

# Results table
df_etas_eid <- data.frame(Database="EID",eta = etas, m.auc=unlist(m.auc),m.thresh=unlist(m.thresh),m.pred=unlist(m.pred),m.hold.out=unlist(m.hold.out),
                 m.auc.all=unlist(m.auc.all),m.thresh.all=unlist(m.thresh.all),m.pred.all=unlist(m.pred.all),m.hold.out.all=unlist(m.hold.out.all))
rownames(df_etas_eid) <- NULL
df_etas_eid

df_etas <- rbind(df_etas_gmpd, df_etas_eid)
saveRDS(df_etas, "non-est-phy-dist-etas-5foldCV.rds")


# Plotting
pdf("eid_phy-dist-eta_5foldCV.pdf")
plot(df_etas_eid$m.auc ~ df_etas_eid$eta, xlab="eta",ylab="Mean AUC across 5 folds", main="AUC for phylogenetic distance transformed by eta")
lines(df_etas_eid$m.auc ~ df_etas_eid$eta)
with(df_etas_eid, abline(v=eta[m.auc==max(m.auc)], col=2, lty=2))
dev.off()





