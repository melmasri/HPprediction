eb.phylo <-
function(phy, heights, a){
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
