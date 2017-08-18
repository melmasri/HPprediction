## Script for the nearest NN

## rm(list= ls())
## Global Variable
SAVE_PARAM = TRUE
SEARCH_NNK = TRUE
ZtZ = FALSE
## DATAFILENAME = 'comEID-PS.RData'
## DATAFILENAME = 'comGMPD.RData'
print(DATAFILENAME)
## source('library.R')
## source('gen.R')
load(DATAFILENAME)
## library(parallel)
library(ape)
library(geiger)
#######################
## Parameters for the independent GGP
## set the correct prior.
## Summary statistics

tree <- read.tree('../mammals.tre')
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(com)])

# Testing all names in hosts in com exist in tree
if(! all(sapply(rownames(com), function(r) r %in% tree$tip.label))) {
    print('Warning! Not all hosts in com exist in tree. Hosts not found in tree will be removed.')
    com <- com[rownames(com)%in%tree$tip.label,]
}

dd = cophenetic(rescale(tree, 'EB', 0))
host.order <- sapply(rownames(dd), function(r) which(r==rownames(com)))
com = com[host.order,]
phy_dist = 1/dd
diag(phy_dist)<-0
com_pa = 1*(com>0)

pairs = cross.validate.fold(com_pa, n=5,1)
tot.gr = length(unique(pairs[,'gr']))


## D = 1*(com>0)
## dist  = D%*%t(D)
## dist  = 1-1/dist
## dist[is.infinite(dist)]<-0
## dist = dist + 1e-4
## diag(dist)<-0
source('../library.R')
source('../gen.R')

if(SEARCH_NNK){
    min.nnk = 2
    max.nnk = nrow(com_pa) - 1
    qartP = floor((max.nnk - min.nnk)/4)
    midP = qartP*2
    qart2P = qartP*3
    pointlist = c(min.nnk, qartP, midP, qart2P, max.nnk)
    auc= c(0,0,0,0,0)
    repeat{
        q = floor((max.nnk - min.nnk)/4)
        pointlist = c(min.nnk, min.nnk + q,  min.nnk + 2*q, min.nnk + 3*q, max.nnk)
        auc1= sapply(pointlist[which(auc==0)], function(nnk)
            mean(unlist(lapply(1:tot.gr ,function(x, pairs, Z, dist, nn.k){
                com_paCross = Z
                com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
                zeros = which(com_paCross==0,arr.ind=TRUE)
                P = com_paCross*0
                Z.in=com_paCross
                if(ZtZ){
                    dist  = Z.in%*%t(Z.in)
                    dist  = 1-1/dist
                    dist[is.infinite(dist)]<-0
                    dist = dist + 1e-4
                    diag(dist)<-0
                }
                nn = apply(dist,1,function(r) order(r, decreasing=TRUE)[1:nn.k])

                z0 = zeros[sample(nrow(zeros)),]
                for(i in 1:nrow(z0))
                    P[z0[i,1],z0[i,2]]<-sum(Z.in[nn[1:nn.k,z0[i,1]],z0[i,2]])
                P = P/nn.k
                P[com_paCross==1]<-1
                roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P, plot=FALSE, bins=400)
                roc$auc
            },pairs=pairs,Z = com_pa, dist=phy_dist,nn.k=nnk))))
        auc[which(auc==0)]<-auc1
        print(rbind(pointlist, auc))
        a =  which.max(auc)
        min.nnk = pointlist[max(1,a-1)]
        max.nnk = pointlist[min(5, a+1)]
        auc[1] = auc[max(1, a-1)]
        auc[5] = auc[min(5,a+1)]
        auc[2]=auc[3]=auc[4]=0
        if (abs(max.nnk  - min.nnk)<4) break
    }    
    nnk = pointlist[a]
    print(sprintf('nnk is: %d', nnk))
}

res = lapply(1:tot.gr ,function(x, pairs, Z, dist, nn.k){
        com_paCross = Z
        com_paCross[pairs[which(pairs[,'gr']==x),c('row', 'col')]]<-0
        zeros = which(com_paCross==0,arr.ind=TRUE)
        P = com_paCross*0
        Z.in=com_paCross
        if(ZtZ){
            dist  = Z.in%*%t(Z.in)
            dist  = 1-1/dist
            dist[is.infinite(dist)]<-0
            dist = dist + 1e-4
            diag(dist)<-0
        }
        nn = apply(dist,1,function(r) order(r, decreasing=TRUE)[1:nn.k])

        z0 = zeros[sample(nrow(zeros)),]
        for(i in 1:nrow(z0))
            P[z0[i,1],z0[i,2]]<-sum(Z.in[nn[1:nn.k,z0[i,1]],z0[i,2]])
        P = P/nn.k
        P[com_paCross==1]<-1

        roc = rocCurves(Z=Z, Z_cross= com_paCross, P=P,
            plot=FALSE, bins=400, all=FALSE)
        tb  = ana.table(Z, com_paCross, roc=roc, plot=FALSE)
        roc.all = rocCurves(Z=Z, Z_cross= com_paCross, P=P,
            plot=FALSE, bins=400, all=TRUE)
        tb.all  = ana.table(Z, com_paCross, roc=roc.all, plot=FALSE)
        
        list(nnk=nn.k, tb = tb, tb.all = tb.all,
             FPR.all = roc.all$roc$FPR, TPR.all=roc.all$roc$TPR,
             FPR = roc$roc$FPR, TPR=roc$roc$TPR)
    },pairs=pairs,Z = com_pa, dist=phy_dist,nn.k=nnk)        


if(SAVE_PARAM)
    save.image(file = 'param.RData')


print(sapply(res, function(r) r$tb))
print(mean(sapply(res, function(r) r$tb$auc)))
##################################################
##################################################


## sapply(c(32:40), function(k){
##     nn = apply(dist,1,function(r) order(r, decreasing=TRUE)[1:k])
##     P = com_paCross*0
##     Z.in=com_paCross
##     z0=which(com_paCross>-1, arr.ind=T)
##     for(i in 1:nrow(z0))
##         P[z0[i,1],z0[i,2]]<-sum(Z.in[nn[1:k,z0[i,1]],z0[i,2]])
##     P = P/k
##     P[com_paCross==1]<-1
##     roc = rocCurves(Z =Z, Z_cross = com_paCross, P=P, plot=FALSE,bins=200)
##     roc$auc
## })

## range(P)
## k=40
## P = com_paCross*0
## Z.in=com_paCross
## nn = apply(dist,1,function(r) order(r, decreasing=TRUE)[1:k])
## P = com_paCross*0
## Z.in=com_paCross
## z0=which(com_paCross>-1, arr.ind=T)
## for(i in 1:nrow(z0))
##     P[z0[i,1],z0[i,2]]<-sum(Z.in[nn[1:k,z0[i,1]],z0[i,2]])
## P = P/k
## roc = rocCurves(Z =Z, Z_cross = com_paCross, P=P, plot=TRUE,bins=200)
## roc$auc
## plot(x=roc$roc$FPR, y=roc$roc$TPR, type='b')
### After teesting GMP k = 33 , EID k= 40
